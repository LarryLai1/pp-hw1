#include <bits/stdc++.h>
#include <fstream>
#include <boost/functional/hash.hpp>
#define MAXSIZE 256
#define mkp make_pair
#define pii pair<int, int>
using namespace std;

/*
o: The player stepping on a regular tile.
O: The player stepping on a target tile.
x: A box on a regular tile.
X: A box on a target tile.
  (space): Nothing on a regular tile.
.: Nothing on a target tile.
#: Wall.
@: A fragile tile where only the player can step on. The boxes are not allowed to be put on it.
!: The player stepping on a fragile tile.
*/

// Directions: up, down, left, right
vector<pair<int, int>> dirs = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};
vector<char> dirChar = {'W','S','A','D'};

vector<string> tmp;
int countRound = 0, countRound_rev = 0;
int moveRound = 0;
int totalr, totalc;
bitset<MAXSIZE> grid = 0, goalPositions = 0, fragPositions = 0;
bitset<MAXSIZE> simpleDeadStateMask = 0;
unordered_map<bitset<MAXSIZE>, bool> isFrozenCache;

struct State {
    pii agentPos; // (row, col)
    bitset<MAXSIZE> cc; // connected components of agentPos
    bitset<MAXSIZE> boxPositions; // bitmask of box positions
    pair<State*, int> path; // (previous state, direction taken to reach this state)

    State() : path({nullptr, 0}) {}
    
    State(pii agent, bitset<MAXSIZE> connComp, bitset<MAXSIZE> boxes, 
          pair<State*, int> last) 
        : agentPos(agent), cc(connComp), boxPositions(boxes), path(last) {}

    bool operator==(const State& other) const {
        return boxPositions == other.boxPositions && agentPos == other.agentPos;
    }
};

// // Comparator for priority queue (min-heap based on fCost)
// struct StateComparator {
//     bool operator()(State* a, State* b) const {
//         if (a->fCost != b->fCost) return a->fCost > b->fCost;
//         return a->hCost > b->hCost; // tie-breaker: prefer lower heuristic
//     }
// };

bitset<MAXSIZE> bitLeft(bitset<MAXSIZE> mp, bitset<MAXSIZE> mask){
    return (mp >> 1) & ~(mask);
}
bitset<MAXSIZE> bitRight(bitset<MAXSIZE> mp, bitset<MAXSIZE> mask){
    return (mp << 1) & ~(mask);
}
bitset<MAXSIZE> bitUp(bitset<MAXSIZE> mp, bitset<MAXSIZE> mask){
    return (mp >> totalc) & ~(mask);
}
bitset<MAXSIZE> bitDown(bitset<MAXSIZE> mp, bitset<MAXSIZE> mask){
    return (mp << totalc) & ~(mask);
}

void printBitset(bitset<MAXSIZE>& bs, const pii& agentPos) {
    for (int i = 0; i < totalr; ++i) {
        for (int j = 0; j < totalc; ++j) {
            int idx = i * totalc + j;
            if (agentPos.first == i && agentPos.second == j) {
                cout << 'A'; // Agent position
            }
            else if (tmp[i][j] == '#') {
                cout << '#'; // Wall
            }
            else {
                cout << (bs[idx] ? '1' : ' ');
            }
        }
        cout << endl;
    }
    cout << endl;
}

// boxPositions, agentPos
pair<bitset<MAXSIZE>, pii> readFileToVector(const string& filename) {
    ifstream file(filename);
    bitset<MAXSIZE> boxPositions = 0;
    pii agentPos;
    
    string line;
    // 逐行讀取檔案內容
    while (getline(file, line)) {
        tmp.push_back(line);
    }
    file.close();
    totalr = tmp.size();
    totalc = tmp[0].size();

    for (int i = 0; i < totalr; ++i) {
        for (int j = 0; j < totalc; ++j) {
            char cell = tmp[i][j];
            int idx = i * totalc + j;
            switch (cell){
                case 'o':
                    agentPos = mkp(i, j);
                    continue;
                case 'O':
                    agentPos = mkp(i, j);
                    goalPositions.set(idx);
                    continue;
                case '@':
                    fragPositions.set(idx);
                    continue;
                case '!':
                    agentPos = mkp(i, j);
                    fragPositions.set(idx);
                    continue;
                case 'x':
                    boxPositions.set(idx);
                    continue;
                case 'X':
                    boxPositions.set(idx);
                    goalPositions.set(idx);
                    continue;
                case '.':
                    goalPositions.set(idx);
                    continue;
                case '#':
                    grid.set(idx);
                    continue;
                default:
                    continue;
            }
        }
    }
    return mkp(boxPositions, agentPos);
}

string agentGoTo(const pii& agentPos, const pii& targetPos, bitset<MAXSIZE> boxPositions) {
    if (agentPos == targetPos) return "";
    int n = totalr;
    int m = totalc;
    int tr = targetPos.first, tc = targetPos.second;
    vector<vector<bool>> visited(n, vector<bool>(m, false));
    queue<tuple<int, int, string>> q;
    q.push({agentPos.first, agentPos.second, ""});
    visited[agentPos.first][agentPos.second] = true;
    while (!q.empty()) {
        auto [r, c, path] = q.front(); q.pop();
        if (r == tr && c == tc) {
            return path;
        }
        for (int d = 0; d < 4; ++d) {
            int nr = r + dirs[d].first;
            int nc = c + dirs[d].second;
            if (nr < 0 || nr >= n || nc < 0 || nc >= m) continue;
            bitset<MAXSIZE> mask = 0;
            mask.set(nr * m + nc);
            // Check if wall or box
            if (((grid & mask) | (boxPositions & mask)).any()) continue;
            if (!visited[nr][nc]) {
                visited[nr][nc] = true;
                q.push({nr, nc, path + dirChar[d]});
            }
        }
    }
    return "N"; // No path found
}

bitset<MAXSIZE> connectedComponent(const pii& startPos, bitset<MAXSIZE> mask) {
    int n = totalr;
    int m = totalc;
    bitset<MAXSIZE> start_mask = 0;
    start_mask.set(startPos.first * m + startPos.second);
    bitset<MAXSIZE> last = 0;
    while (last != start_mask) {
        last = start_mask;
        start_mask |= bitLeft(start_mask, mask);
        start_mask |= bitRight(start_mask, mask);
        start_mask |= bitUp(start_mask, mask);
        start_mask |= bitDown(start_mask, mask);
    }
    return start_mask;
}

// 0: reachable, 1: simple dead state
void simpleDeadlockList(){
    queue<pair<int, int>> goalQueue;

    vector<vector<int>> simpleDeadState(tmp.size(), vector<int>(tmp[0].size(), 0));
    for (int i = 0; i < totalr; ++i) {
        for (int j = 0; j < totalc; ++j) {
            char cell = tmp[i][j];
            if (cell == '#') continue; // wall
            if (cell == '.' || cell == 'O' || cell == 'X') {
                goalQueue.push(mkp(i, j));
                continue; // target
            }
            // empty space
            simpleDeadState[i][j] = 1;
        }
    }

    while (goalQueue.size() > 0) {
        pair<int, int> goal = goalQueue.front(); goalQueue.pop();
        int r = goal.first;
        int c = goal.second;
        // Check all four directions
        for (const auto& dir : dirs) {
            int nr = r + dir.first;
            int nc = c + dir.second;
            if (simpleDeadState[nr][nc] == 1 && tmp[nr+dir.first][nc+dir.second] != '#') {
                simpleDeadState[nr][nc] = 0; // Mark as reachable
                goalQueue.push(mkp(nr, nc));
            }
        }
    }

    for (int i = 0; i < totalr; ++i) {
        for (int j = 0; j < totalc; ++j) {
            if (simpleDeadState[i][j] == 1) {
                simpleDeadStateMask.set(i * totalc + j);
            }
        }
    }
}

// return up, down, left, right possible moves
tuple<bitset<MAXSIZE>, bitset<MAXSIZE>, bitset<MAXSIZE>
    , bitset<MAXSIZE>> getPossibleMoves(const bitset<MAXSIZE>& boxPositions) {
    auto leftMoves = bitRight(bitLeft(boxPositions, grid | fragPositions | boxPositions | simpleDeadStateMask), 0) &
            ~(bitLeft(bitRight(boxPositions, ~(grid | boxPositions)), 0));
    auto rightMoves = bitLeft(bitRight(boxPositions, grid | fragPositions | boxPositions | simpleDeadStateMask), 0) &
            ~(bitRight(bitLeft(boxPositions, ~(grid | boxPositions)), 0));
    auto upMoves = bitDown(bitUp(boxPositions, grid | fragPositions | boxPositions | simpleDeadStateMask), 0) &
            ~(bitUp(bitDown(boxPositions, ~(grid | boxPositions)), 0));
    auto downMoves = bitUp(bitDown(boxPositions, grid | fragPositions | boxPositions | simpleDeadStateMask), 0) &
            ~(bitDown(bitUp(boxPositions, ~(grid | boxPositions)), 0));
    return {upMoves, downMoves, leftMoves, rightMoves};
}

bool isDeadFrozen(const bitset<MAXSIZE>& boxPositions){
    if (isFrozenCache.find(boxPositions) != isFrozenCache.end()) {
        return isFrozenCache[boxPositions];
    }
    bitset<MAXSIZE> leftMoves = 0, rightMoves = 0, upMoves = 0, downMoves = 0;
    bitset<MAXSIZE> anyMove = 0, tmpbox = boxPositions;
    while (true){
        leftMoves = bitRight(bitLeft(tmpbox, grid | fragPositions | tmpbox | simpleDeadStateMask), 0) &
                ~(bitLeft(bitRight(tmpbox, ~(grid | tmpbox)), 0));
        rightMoves = bitLeft(bitRight(tmpbox, grid | fragPositions | tmpbox | simpleDeadStateMask), 0) &
                ~(bitRight(bitLeft(tmpbox, ~(grid | tmpbox)), 0));
        upMoves = bitDown(bitUp(tmpbox, grid | fragPositions | tmpbox | simpleDeadStateMask), 0) &
                ~(bitUp(bitDown(tmpbox, ~(grid | tmpbox)), 0));
        downMoves = bitUp(bitDown(tmpbox, grid | fragPositions | tmpbox | simpleDeadStateMask), 0) &
                ~(bitDown(bitUp(tmpbox, ~(grid | tmpbox)), 0));
        auto newAnyMove = leftMoves | rightMoves | upMoves | downMoves;
        if (newAnyMove == anyMove) break;
        anyMove = newAnyMove;
        tmpbox &= ~anyMove;
    }
    auto deadFrozen = (tmpbox & ~goalPositions).any();
    isFrozenCache[boxPositions] = deadFrozen;
    return deadFrozen;
}

// locate movable boxes of all directions considering connected component of agent
// 0: up, 1: down, 2: left, 3: right
vector<bitset<MAXSIZE>> getRealMoves(const bitset<MAXSIZE>& boxPos, 
    const bitset<MAXSIZE>& cc, tuple<bitset<MAXSIZE>, bitset<MAXSIZE>, bitset<MAXSIZE>, bitset<MAXSIZE>> moves) {
    auto [upMoves, downMoves, leftMoves, rightMoves] = moves;
    bitset<MAXSIZE> left = bitLeft(bitRight(leftMoves, ~cc), 0);
    bitset<MAXSIZE> right = bitRight(bitLeft(rightMoves, ~cc), 0);
    bitset<MAXSIZE> up = bitUp(bitDown(upMoves, ~cc), 0);
    bitset<MAXSIZE> down = bitDown(bitUp(downMoves, ~cc), 0);
    return {up, down, left, right};
}

// return up, down, left, right possible moves
tuple<bitset<MAXSIZE>, bitset<MAXSIZE>, bitset<MAXSIZE>
    , bitset<MAXSIZE>> getPossibleMovesRev(const bitset<MAXSIZE>& boxPositions) {
    auto leftMoves = bitLeft(bitLeft(bitRight(bitRight(boxPositions, grid | fragPositions | boxPositions | simpleDeadStateMask), 
        grid | boxPositions), 0), 0);
    auto rightMoves = bitRight(bitRight(bitLeft(bitLeft(boxPositions, grid | fragPositions | boxPositions | simpleDeadStateMask), 
        grid | boxPositions), 0), 0);
    auto upMoves = bitUp(bitUp(bitDown(bitDown(boxPositions, grid | fragPositions | boxPositions | simpleDeadStateMask), 
        grid | boxPositions), 0), 0);
    auto downMoves = bitDown(bitDown(bitUp(bitUp(boxPositions, grid | fragPositions | boxPositions | simpleDeadStateMask), 
        grid | boxPositions), 0), 0);
    
    return {upMoves, downMoves, leftMoves, rightMoves};
}

bool isDeadFrozenRev(const bitset<MAXSIZE>& boxPositions){
    if (isFrozenCache.find(boxPositions) != isFrozenCache.end()) {
        return isFrozenCache[boxPositions];
    }
    bitset<MAXSIZE> leftMoves = 0, rightMoves = 0, upMoves = 0, downMoves = 0;
    bitset<MAXSIZE> anyMove = 0, tmpbox = boxPositions;
    while (true){
        leftMoves = bitLeft(bitLeft(bitRight(bitRight(tmpbox, grid | fragPositions | tmpbox | simpleDeadStateMask), 
            grid | tmpbox), 0), 0);
        rightMoves = bitRight(bitRight(bitLeft(bitLeft(tmpbox, grid | fragPositions | tmpbox | simpleDeadStateMask), 
            grid | tmpbox), 0), 0);
        upMoves = bitUp(bitUp(bitDown(bitDown(tmpbox, grid | fragPositions | tmpbox | simpleDeadStateMask), 
            grid | tmpbox), 0), 0);
        downMoves = bitDown(bitDown(bitUp(bitUp(tmpbox, grid | fragPositions | tmpbox | simpleDeadStateMask), 
            grid | tmpbox), 0), 0);
        auto newAnyMove = leftMoves | rightMoves | upMoves | downMoves;
        if (newAnyMove == anyMove) break;
        anyMove = newAnyMove;
        tmpbox &= ~anyMove;
    }
    auto deadFrozen = (tmpbox & ~goalPositions).any();
    isFrozenCache[boxPositions] = deadFrozen;
    return deadFrozen;
}
vector<bitset<MAXSIZE>> getRealMovesRev(const bitset<MAXSIZE>& boxPos, 
    const bitset<MAXSIZE>& cc, tuple<bitset<MAXSIZE>, bitset<MAXSIZE>, bitset<MAXSIZE>, bitset<MAXSIZE>> moves) {
    auto [upMoves, downMoves, leftMoves, rightMoves] = moves;
    bitset<MAXSIZE> left = bitLeft(bitLeft(bitRight(bitRight(leftMoves, 0), 
        ~cc), 0), 0);
    bitset<MAXSIZE> right = bitRight(bitRight(bitLeft(bitLeft(rightMoves, 0), 
        ~cc), 0), 0);
    bitset<MAXSIZE> up = bitUp(bitUp(bitDown(bitDown(upMoves, 0), 
        ~cc), 0), 0);
    bitset<MAXSIZE> down = bitDown(bitDown(bitUp(bitUp(downMoves, 0), 
        ~cc), 0), 0);
    return {up, down, left, right};
}


bool isIn(const State& state, unordered_map<bitset<MAXSIZE>, 
    pair<bitset<MAXSIZE>, vector<State*>>> &visited) {
    auto found = visited.find(state.boxPositions);
    if (found == visited.end()) return false;
    bitset<MAXSIZE> cc = found->second.first;
    return cc[state.agentPos.first * totalc + state.agentPos.second];
}

State* getState(State* state, unordered_map<bitset<MAXSIZE>, 
    pair<bitset<MAXSIZE>, vector<State*>>> &visited) {
    auto found = visited.find(state->boxPositions);
    auto vec = visited[state->boxPositions].second;
    for (auto s : vec) {
        if (s->cc[state->agentPos.first * totalc + state->agentPos.second]) {
            return s;
        }
    }
    return nullptr;
}

void push_in(State* state, vector<State*> &nv, 
    unordered_map<bitset<MAXSIZE>, pair<bitset<MAXSIZE>, vector<State*>>> &visited) {
    visited[state->boxPositions].first |= state->cc;
    visited[state->boxPositions].second.push_back(state);
    nv.push_back(state);
}

string recoverPath(State* initialState, State* midState, State* midState_rev) {
    string result;
    vector<pair<pii, int>> steps; // (box position, direction)
    State* curstate = midState;
    while (curstate != nullptr && curstate->path.first != nullptr) {
        steps.push_back({curstate->agentPos, curstate->path.second});
        curstate = (curstate->path).first;
    }

    reverse(steps.begin(), steps.end());
    pii currentPos = initialState->agentPos;
    auto currentBoxPos = initialState->boxPositions;
    for (const auto& step : steps) {
        auto [nextPos, dir] = step;
        pii targetPos = {nextPos.first - dirs[dir].first, nextPos.second - dirs[dir].second}; // Position behind the box
        string pathToTarget = agentGoTo(currentPos, targetPos, currentBoxPos);
        if (pathToTarget == "N") {
            cout << "Something went wrong, no path found!" << endl;
            return "N"; // No valid path found
        }
        result += pathToTarget; // Move agent to the position behind the box
        result += dirChar[dir]; // Push the box
        currentPos = nextPos;
        // Update currentBoxPos after each push
        currentBoxPos.reset(nextPos.first * totalc + nextPos.second);
        currentBoxPos.set((nextPos.first + dirs[dir].first) * totalc + nextPos.second + dirs[dir].second);
    }
    curstate = midState_rev;
    while (curstate != nullptr && curstate->path.first != nullptr) {
        result += agentGoTo(currentPos, curstate->agentPos, curstate->boxPositions);
        // cout << result << endl;
        result += dirChar[(curstate->path.second)]; // Reverse direction
        currentPos = {curstate->agentPos.first + dirs[curstate->path.second].first, 
                      curstate->agentPos.second + dirs[curstate->path.second].second};
        curstate = (curstate->path).first;
        // cout << "CurPos: (" << currentPos.first << ", " << currentPos.second << ")\n";
    }
    return result;
}

void bfs(State* initialState, State*endState, vector<State*> &v, vector<State*> &nv, 
    unordered_map<bitset<MAXSIZE>, pair<bitset<MAXSIZE>, vector<State*>>> &visited,
    unordered_map<bitset<MAXSIZE>, pair<bitset<MAXSIZE>, vector<State*>>> &visited_another) {
    for (auto state : v) {
        auto [agentPos, cc, boxPositions, path] = *state;
        auto movableBoxes = getPossibleMoves(boxPositions);
        auto boxMoves = getRealMoves(boxPositions, cc, movableBoxes);
        
        for (int k=0; k<4; ++k){
            bitset<MAXSIZE> curmove = boxMoves[k];
            size_t idx = curmove._Find_first();
            while (idx < curmove.size()) {
                int boxIdx = idx;
                int boxRow = boxIdx / totalc;
                int boxCol = boxIdx % totalc;

                // Compute new box position after move
                int dr = dirs[k].first;
                int dc = dirs[k].second;
                int newRow = boxRow + dr;
                int newCol = boxCol + dc;
                if (newRow < 0 || newRow >= totalr || newCol < 0 || newCol >= totalc) {
                    idx = curmove._Find_next(idx);
                    continue;
                }
                int newBoxIdx = newRow * totalc + newCol;

                // Update box positions
                bitset<MAXSIZE> newBoxPositions = boxPositions;
                newBoxPositions.reset(boxIdx);
                newBoxPositions.set(newBoxIdx);

                // Check for deadlock
                auto deadFrozen = isDeadFrozen(newBoxPositions);
                if (deadFrozen) {
                    idx = curmove._Find_next(idx);
                    continue;
                }

                // Compute agent position (where agent must be to push box)
                pii newAgentPos = mkp(boxRow, boxCol);
                
                // check if in available connected component
                if (isIn({newAgentPos, 0, newBoxPositions, {}}, visited)) {
                    idx = curmove._Find_next(idx);
                    continue;
                }

                // Update connected component for new agent position
                auto newCC = connectedComponent(newAgentPos, grid | newBoxPositions);
                
                // Create new state
                State* newState = new State(newAgentPos, newCC, newBoxPositions, {state, k});
                
                // Check for solution
                if (isIn(*newState, visited_another)) {
                    State *foundState = getState(newState, visited_another);
                    cout << recoverPath(initialState, newState, foundState);
                    cout << endl;
                    exit(0);
                }

                push_in(newState, nv, visited);
                idx = curmove._Find_next(idx);
            }
        }
        countRound++;
    }
}

void bfs_rev(State* initialState, State*endState, vector<State*> &v, vector<State*> &nv, 
    unordered_map<bitset<MAXSIZE>, pair<bitset<MAXSIZE>, vector<State*>>> &visited, 
    unordered_map<bitset<MAXSIZE>, pair<bitset<MAXSIZE>, vector<State*>>> &visited_another) {
    for (auto state : v) {
        auto [agentPos, cc, boxPositions, path] = *state;
        auto movableBoxes = getPossibleMovesRev(boxPositions);
        auto boxMoves = getRealMovesRev(boxPositions, cc, movableBoxes);
        
        for (int k=0; k<4; ++k){
            bitset<MAXSIZE> curmove = boxMoves[k];
            size_t idx = curmove._Find_first();
            while (idx < curmove.size()) {
                int boxIdx = idx;
                int boxRow = boxIdx / totalc;
                int boxCol = boxIdx % totalc;

                // Compute new box position after move
                int dr = dirs[k].first;
                int dc = dirs[k].second;
                int newRow = boxRow - dr;
                int newCol = boxCol - dc;
                if (newRow < 0 || newRow >= totalr || newCol < 0 || newCol >= totalc) {
                    idx = curmove._Find_next(idx);
                    continue;
                }
                int newBoxIdx = newRow * totalc + newCol;

                // Update box positions
                bitset<MAXSIZE> newBoxPositions = boxPositions;
                newBoxPositions.reset(boxIdx);
                newBoxPositions.set(newBoxIdx);

                // Check for deadlock
                auto deadFrozen = isDeadFrozenRev(newBoxPositions);
                if (deadFrozen) {
                    idx = curmove._Find_next(idx);
                    continue;
                }

                // Compute agent position (where agent must be to push box)
                int agentRow = boxRow - 2*dr;
                int agentCol = boxCol - 2*dc;
                pii newAgentPos = mkp(agentRow, agentCol);
                
                // check if in available connected component
                if (isIn({newAgentPos, 0, newBoxPositions, {}}, visited)) {
                    idx = curmove._Find_next(idx);
                    continue;
                }

                // Update connected component for new agent position
                auto newCC = connectedComponent(newAgentPos, grid | newBoxPositions);
                
                // Create new state
                State* newState = new State(newAgentPos, newCC, newBoxPositions, {state, k});

                // Check for solution
                if (isIn(*newState, visited_another)) {
                    State *foundState = getState(newState, visited_another);
                    cout << recoverPath(initialState, foundState, newState);
                    cout << endl;
                    exit(0);
                }

                push_in(newState, nv, visited);
                idx = curmove._Find_next(idx);
            }
        }
        countRound_rev++;
    }
}

int main(int argc, char* argv[]) {
    // 檢查命令列參數
    if (argc != 2) {
        cerr << "使用方式: " << argv[0] << " <檔案路徑>" << endl;
        cerr << "範例: ./worker 01.txt" << endl;
        return 1;
    }
    // read file
    string filename = argv[1];
    auto [boxPositions, agentPos] = readFileToVector(filename);
    // cout << "totalr: " << totalr << ", totalc: " << totalc << endl;
    // cout << "==========================" << endl;
    auto cc = connectedComponent(agentPos, grid | boxPositions);
    simpleDeadlockList();

    State* initialState = new State(agentPos, cc, boxPositions, {nullptr, 0});
    State* endState = new State({-1, -1}, ~(bitset<MAXSIZE>)0, goalPositions, {nullptr, 0});
    vector<State*> v, nv;
    unordered_map<bitset<MAXSIZE>, pair<bitset<MAXSIZE>, vector<State*>>> visited;
    vector<State*> v_rev, nv_rev;
    unordered_map<bitset<MAXSIZE>, pair<bitset<MAXSIZE>, vector<State*>>> visited_rev;

    push_in(initialState, v, visited);
    push_in(endState, v_rev, visited_rev);
    #ifdef DEBUG
        printBitset(simpleDeadStateMask, {-1, -1});
    #endif

    while (!v.empty()){
        #ifdef DEBUG
        cout << "Round " << moveRound << ", Forward frontier size: " << v.size() << ", Reverse frontier size: " << v_rev.size()
             << ", Forward count: "<< countRound << ", Reverse count: " << countRound_rev << endl;
        #endif
        bfs(initialState, endState, v, nv, visited, visited_rev);
        v = nv;
        nv.clear();
        bfs_rev(initialState, endState, v_rev, nv_rev, visited_rev, visited);
        v_rev = nv_rev;
        nv_rev.clear();
        moveRound++;
    }
    cout << "No solution found." << endl;

    return 0;
}