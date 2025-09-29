#include <bits/stdc++.h>
#include <fstream>
#include <assert.h>
#include <boost/functional/hash.hpp>
#include <signal.h>
#include <unistd.h>
#define MAXSIZE 256

void alarm_handler(int signum) {
    std::cerr << "Time limit exceeded." << std::endl;
    exit(1);
}

struct AlarmSetter {
    AlarmSetter() {
        signal(SIGALRM, alarm_handler);
        alarm(180);
    }
} alarmSetterInstance;
#define mkp std::make_pair
#define pii std::pair<int, int>
#define deadlocklimit 10000

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
std::vector<std::pair<int, int>> dirs = {{-1, 0}, {1, 0}, {0, -1}, {0, 1}};
std::vector<char> dirChar = {'W','S','A','D'};
struct State {
    pii agentPos; // (row, col)
    std::bitset<MAXSIZE> cc; // connected components of agentPos
    std::bitset<MAXSIZE> boxPositions; // bitmask of box positions
    std::vector<std::pair<pii, int>> path; // path taken to reach this state

    bool operator==(const State& other) const {
        return boxPositions == other.boxPositions && agentPos == other.agentPos;
    }
};

std::vector<std::string> tmp;
int count = 0;
int moveRound = 0;
int totalr, totalc;
std::bitset<MAXSIZE> grid = 0, goalPositions = 0, fragPositions = 0;
std::bitset<MAXSIZE> simpleDeadStateMask = 0;
std::unordered_map<std::bitset<MAXSIZE>, bool> isFrozenCache;

// TODO: it may not be correct
std::bitset<MAXSIZE> bitLeft(std::bitset<MAXSIZE> mp, std::bitset<MAXSIZE> mask){
    return (mp >> 1) & ~(mask);
}
std::bitset<MAXSIZE> bitRight(std::bitset<MAXSIZE> mp, std::bitset<MAXSIZE> mask){
    return (mp << 1) & ~(mask);
}
std::bitset<MAXSIZE> bitUp(std::bitset<MAXSIZE> mp, std::bitset<MAXSIZE> mask){
    return (mp >> totalc) & ~(mask);
}
std::bitset<MAXSIZE> bitDown(std::bitset<MAXSIZE> mp, std::bitset<MAXSIZE> mask){
    return (mp << totalc) & ~(mask);
}


// boxPositions, agentPos
std::pair<std::bitset<MAXSIZE>, pii> readFileToVector(const std::string& filename) {
    std::ifstream file(filename);
    std::bitset<MAXSIZE> boxPositions = 0;
    pii agentPos;

    assert (file.is_open());
    
    std::string line;
    // 逐行讀取檔案內容
    while (std::getline(file, line)) {
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

std::string agentGoTo(const pii& agentPos, const pii& targetPos, std::bitset<MAXSIZE> boxPositions) {
    if (agentPos == targetPos) return "";
    int n = totalr;
    int m = totalc;
    int tr = targetPos.first, tc = targetPos.second;
    std::vector<std::vector<bool>> visited(n, std::vector<bool>(m, false));
    std::queue<std::tuple<int, int, std::string>> q;
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
            std::bitset<MAXSIZE> mask = 0;
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

std::bitset<MAXSIZE> connectedComponent(const pii& startPos, std::bitset<MAXSIZE> mask) {
    int n = totalr;
    int m = totalc;
    std::bitset<MAXSIZE> start_mask = 0;
    start_mask.set(startPos.first * m + startPos.second);
    std::bitset<MAXSIZE> last = 0;
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
    std::queue<std::pair<int, int>> goalQueue;

    std::vector<std::vector<int>> simpleDeadState(tmp.size(), std::vector<int>(tmp[0].size(), 0));
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
        std::pair<int, int> goal = goalQueue.front(); goalQueue.pop();
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
std::tuple<std::bitset<MAXSIZE>, std::bitset<MAXSIZE>, std::bitset<MAXSIZE>
    , std::bitset<MAXSIZE>> getPossibleMoves(const std::bitset<MAXSIZE>& boxPositions) {
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

bool isDeadFrozen(const std::bitset<MAXSIZE>& boxPositions){
    if (isFrozenCache.find(boxPositions) != isFrozenCache.end()) {
        return isFrozenCache[boxPositions];
    }
    std::bitset<MAXSIZE> leftMoves = 0, rightMoves = 0, upMoves = 0, downMoves = 0;
    std::bitset<MAXSIZE> anyMove = 0, tmpbox = boxPositions;
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
// TODO: this may be wrong
// 0: up, 1: down, 2: left, 3: right
std::vector<std::bitset<MAXSIZE>> getRealMoves(const std::bitset<MAXSIZE>& boxPos, 
    const std::bitset<MAXSIZE>& cc, std::tuple<std::bitset<MAXSIZE>, std::bitset<MAXSIZE>, std::bitset<MAXSIZE>, std::bitset<MAXSIZE>> moves) {
    auto [upMoves, downMoves, leftMoves, rightMoves] = moves;
    std::bitset<MAXSIZE> left = bitLeft(bitRight(leftMoves, ~cc), 0);
    std::bitset<MAXSIZE> right = bitRight(bitLeft(rightMoves, ~cc), 0);
    std::bitset<MAXSIZE> up = bitUp(bitDown(upMoves, ~cc), 0);
    std::bitset<MAXSIZE> down = bitDown(bitUp(downMoves, ~cc), 0);
    return {up, down, left, right};
}

bool isIn(const State& state, std::unordered_map<std::bitset<MAXSIZE>, std::bitset<MAXSIZE>>& visited) {
    std::bitset<MAXSIZE> cc = visited[state.boxPositions];
    return cc[state.agentPos.first * totalc + state.agentPos.second];
}

void push_in(const State& state, std::vector<State> &v, std::unordered_map<std::bitset<MAXSIZE>, 
    std::bitset<MAXSIZE>> &visited) {
    visited[state.boxPositions] |= state.cc;
    v.push_back(state);
}

void printBitset(const std::bitset<MAXSIZE>& bs, const pii& agentPos) {
    for (int i = 0; i < totalr; ++i) {
        for (int j = 0; j < totalc; ++j) {
            int idx = i * totalc + j;
            if (agentPos.first == i && agentPos.second == j) {
                std::cout << 'A'; // Agent position
            }
            else if (tmp[i][j] == '#') {
                std::cout << '#'; // Wall
            }
            else {
                std::cout << (bs[idx] ? '1' : ' ');
            }
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

std::string recoverPath(const State& initialState, const State& finalState) {
    std::string result;
    pii currentPos = initialState.agentPos;
    auto currentBoxPos = initialState.boxPositions;
    for (const auto& step : finalState.path) {
        auto [targetPos, dir] = step;
        std::string pathToTarget = agentGoTo(currentPos, targetPos, currentBoxPos);
        if (pathToTarget == "N") {
            std::cout << "Something went wrong, no path found!" << std::endl;
            return "N"; // No valid path found
        }
        result += pathToTarget; // Move agent to the position behind the box
        result += dirChar[dir]; // Push the box
        currentPos = {targetPos.first + dirs[dir].first, targetPos.second + dirs[dir].second}; // Update agent position to the box's new position
        // Update currentBoxPos after each push
        currentBoxPos.reset(currentPos.first * totalc + currentPos.second);
        currentBoxPos.set((currentPos.first + dirs[dir].first) * totalc + currentPos.second + dirs[dir].second);
    }
    return result;
}

void bfs(const State& initialState) {
    std::vector<State> v, nv;
    // Initialize with movable boxes for initial state
    std::unordered_map<std::bitset<MAXSIZE>, std::bitset<MAXSIZE>> visited;

    push_in(initialState, v, visited);
    
    while (!v.empty()) {
        #ifdef DEBUG
            std::cout << "Round: " << moveRound << ", States: " << v.size() << ", Total: " << count << std::endl;
        #endif
        for (const auto& state: v){
            auto [agentPos, cc, boxPositions, path] = state;
            auto movableBoxes = getPossibleMoves(boxPositions);
            auto boxMoves = getRealMoves(boxPositions, state.cc, movableBoxes);
            for (int k=0; k<4; ++k){
                std::bitset<MAXSIZE> curmove = boxMoves[k];
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
                    std::bitset<MAXSIZE> newBoxPositions = boxPositions;
                    newBoxPositions.reset(boxIdx);
                    newBoxPositions.set(newBoxIdx);

                    // Check for deadlock
                    auto deadFrozen = isDeadFrozen(newBoxPositions);
                    if (deadFrozen) {
                        idx = curmove._Find_next(idx);
                        continue;
                    }

                    // Compute agent position (where agent must be to push box)
                    int agentRow = boxRow - dr;
                    int agentCol = boxCol - dc;
                    pii newAgentPos = mkp(boxRow, boxCol);
                    pii requiredAgentPos = mkp(agentRow, agentCol);

                    // Update path
                    std::vector<std::pair<pii, int>> newPath = state.path;
                    newPath.push_back({requiredAgentPos, k});
                    
                    // check if in available connected component
                    if (isIn({newAgentPos, 0, newBoxPositions, {}}, visited)) {
                        idx = curmove._Find_next(idx);
                        continue;
                    }
                    // Update connected component for new agent position
                    auto newCC = connectedComponent(newAgentPos, grid | newBoxPositions);
                    State newState = {newAgentPos, newCC, newBoxPositions, newPath};

                    // Check for solution
                    if ((newBoxPositions & ~goalPositions).none()) {
                        std::cout << recoverPath(initialState, newState);
                        std::cout << std::endl;
                        exit(0);
                    }
                    push_in(newState, nv, visited);
                    idx = curmove._Find_next(idx);
                }
            }
            count++;
        }
        v = nv;
        nv.clear();
        moveRound++;
    }
}

int main(int argc, char* argv[]) {
    // 檢查命令列參數
    if (argc != 2) {
        std::cerr << "使用方式: " << argv[0] << " <檔案路徑>" << std::endl;
        std::cerr << "範例: ./worker 01.txt" << std::endl;
        return 1;
    }
    // read file
    std::string filename = argv[1];
    auto [boxPositions, agentPos] = readFileToVector(filename);
    // std::cout << "totalr: " << totalr << ", totalc: " << totalc << std::endl;
    // std::cout << "==========================" << std::endl;
    auto cc = connectedComponent(agentPos, grid | boxPositions);
    simpleDeadlockList();
    bfs({agentPos, cc, boxPositions, {}});

    return 0;
}