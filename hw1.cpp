#include <bits/stdc++.h>
#include <fstream>
#include <assert.h>
#include <boost/functional/hash.hpp>
#include <signal.h>
#include <unistd.h>
#define MAXSIZE 100

/*
To fix:
save movable boxes in State, StatePos in std::set<pii> instead of std::vector<pii>
also change the output type of getMovableBoxes.first to std::set<pii>
change isIn, pushin
do getMovableBoxes and check deadfrozen before pushing in (need to draw curgrid everytime)
that is to say, only need to do getAllBoxMoves after popping from v
*/

void alarm_handler(int signum) {
    std::cerr << "Time limit exceeded." << std::endl;
    exit(1);
}

struct AlarmSetter {
    AlarmSetter() {
        signal(SIGALRM, alarm_handler);
        alarm(60);
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

struct StatePos{
    std::bitset<MAXSIZE> cc; // all connected components of agentPos
    std::bitset<MAXSIZE> boxPositions; // (row, col)
    // Hash function for StatePos using only boxPositions
    bool operator==(const StatePos& other) const {
        return boxPositions == other.boxPositions;
    }
};

namespace std {
    template <>
    struct hash<StatePos> {
        size_t operator()(const StatePos& s) const {
            // Only hash boxPositions (std::bitset<MAXSIZE>)
            return std::hash<std::bitset<MAXSIZE>>()(s.boxPositions);
        }
    };
}

// Cache for storing movableBoxes under different boxPositions
std::vector<std::string> tmp;
int count = 0;
int moveRound = 0;
int totalr, totalc;
std::bitset<MAXSIZE> grid = 0, goalPositions = 0, fragPositions = 0;
std::bitset<MAXSIZE> simpleDeadStateMask = 0;

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
    // std::cout << "doing bitDown" << std::endl;
    // std::cout << "mp: \t\t\t\t" << mp << std::endl;
    // std::cout << "mask: \t\t\t\t" << mask << std::endl;
    // std::cout << "(mp << totalc): \t\t" << (mp << totalc) << std::endl;
    // std::cout << "(mp << totalc) & ~(mask): \t" << ((mp << totalc) & ~(mask)) << std::endl;
    // std::cout << "-----------------------------------------" << std::endl;
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
            std::bitset<MAXSIZE> mask = 1u << (nr * m + nc);
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

// bool hasPath(const pii& agentPos, const pii& targetPos, std::bitset<MAXSIZE> mask) {
//     int n = totalr;
//     int m = totalc;
//     if (agentPos == targetPos) return true;
//     std::bitset<MAXSIZE> target_mask = 1u << (targetPos.first * m + targetPos.second);
//     std::bitset<MAXSIZE> start_mask = 1u << (agentPos.first * m + agentPos.second);
//     std::bitset<MAXSIZE> last = 0;
//     while (last != start_mask) {
//         last = start_mask;
//         start_mask |= bitLeft(start_mask, mask);
//         start_mask |= bitRight(start_mask, mask);
//         start_mask |= bitUp(start_mask, mask);
//         start_mask |= bitDown(start_mask, mask);
//         if ((start_mask & target_mask).any()) return true; // Path found
//     }
//     return false; // No path found
// }

std::bitset<MAXSIZE> connectedComponent(const pii& startPos, std::bitset<MAXSIZE> mask) {
    int n = totalr;
    int m = totalc;
    std::bitset<MAXSIZE> start_mask = 1u << (startPos.first * m + startPos.second);
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
std::pair<std::tuple<std::bitset<MAXSIZE>, std::bitset<MAXSIZE>, std::bitset<MAXSIZE>, std::bitset<MAXSIZE>>, 
    bool> getPossibleMoves(const std::bitset<MAXSIZE>& boxPositions) {
    auto leftMoves = bitRight(bitLeft(boxPositions, grid | fragPositions | boxPositions | simpleDeadStateMask), 0);
    auto rightMoves = bitLeft(bitRight(boxPositions, grid | fragPositions | boxPositions | simpleDeadStateMask), 0);
    auto upMoves = bitDown(bitUp(boxPositions, grid | fragPositions | boxPositions | simpleDeadStateMask), 0);
    auto downMoves = bitUp(bitDown(boxPositions, grid | fragPositions | boxPositions | simpleDeadStateMask), 0);
    bool isDeadlocked = ((leftMoves | rightMoves | upMoves | downMoves) == 0) && (boxPositions & ~goalPositions).any();
    return {std::make_tuple(upMoves, downMoves, leftMoves, rightMoves), isDeadlocked};
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

bool isIn(const State& state, const std::unordered_set<StatePos> &visited){
    StatePos statePos;
    statePos.boxPositions = state.boxPositions;
    auto found = visited.find(statePos);
    if (found != visited.end()) {
        auto cc = found->cc;
        if (cc[state.agentPos.first * totalc + state.agentPos.second]) {
            return true;
        }
    }
    return false;
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
                std::cout << (bs[idx] ? '1' : '0');
            }
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void pushin(const State& state, std::vector<State> &v, std::unordered_set<StatePos> &visited) {
    #ifdef DEBUG
        std::cout << "Pushing in new state" << std::endl;
        std::cout << "Agent position: (" << state.agentPos.first << ", " << state.agentPos.second << ")" << std::endl;
        std::cout << "Box positions: \n";
        printBitset(state.boxPositions, state.agentPos);
    #endif
    StatePos statePos;
    statePos.cc = state.cc;
    statePos.boxPositions = state.boxPositions;
    auto found = visited.find(statePos);
    if (found != visited.end()) {
        auto cc = found->cc;
        if (cc[state.agentPos.first * totalc + state.agentPos.second]) {
            return;
        }
        statePos.cc |= found->cc;
        v.push_back(state);
        visited.insert(statePos);
        visited.erase(found);
        return;
    }
    // not found, insert new state
    visited.insert(statePos);
    v.push_back(state);
}

std::string recoverPath(const State& initialState, const State& finalState) {
    std::string result;
    pii currentPos = initialState.agentPos;
    auto currentBoxPos = initialState.boxPositions;
    // std::cout << "Recovering path..." << std::endl;
    // for (const auto& step : finalState.path) {
        // auto [targetPos, dir] = step;
        // std::cout << "Step to (" << targetPos.first << ", " << targetPos.second << ") with push " << dir << std::endl;
    // }
    // std::cout << "-------------------------" << std::endl;
    for (const auto& step : finalState.path) {
        auto [targetPos, dir] = step;
        // std::cout << "Current agent position: (" << currentPos.first << ", " << currentPos.second << ")" << std::endl;
        // std::cout << result << std::endl;
        std::string pathToTarget = agentGoTo(currentPos, targetPos, currentBoxPos);
        if (pathToTarget == "N") {
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
    std::unordered_set<StatePos> visited;

    pushin(initialState, v, visited);
    
    while (!v.empty()) {
        std::cout << "==========================" << std::endl;
        std::cout << "Round: " << moveRound << ", States: " << v.size() << ", Total: " << count << std::endl;
        for (const auto& state: v){
            auto [agentPos, cc, boxPositions, path] = state;
            auto [movableBoxes, deadFrozen] = getPossibleMoves(boxPositions);
            #ifdef DEBUG
                std::cout << "Box positions: \t\t\t" << boxPositions << std::endl;
                std::cout << "Connected component: \t\t" << cc << std::endl;
                std::cout << "Up Possible Moves: \t\t" << std::get<0>(movableBoxes) << std::endl;
                std::cout << "Down Possible Moves: \t\t" << std::get<1>(movableBoxes) << std::endl;
                std::cout << "Left Possible Moves: \t\t" << std::get<2>(movableBoxes) << std::endl;
                std::cout << "Right Possible Moves: \t\t" << std::get<3>(movableBoxes) << std::endl;
            #endif
            if (deadFrozen) continue;
            #ifdef DEBUG
                std::cout << "Current state with agent at (" << agentPos.first << ", " << agentPos.second << ")" << std::endl;
            #endif
            auto boxMoves = getRealMoves(boxPositions, state.cc, movableBoxes);
            #ifdef DEBUG
                std::cout << "Up Real Moves: \t\t" << boxMoves[0] << std::endl;
                std::cout << "Down Real Moves: \t" << boxMoves[1] << std::endl;
                std::cout << "Left Real Moves: \t" << boxMoves[2] << std::endl;
                std::cout << "Right Real Moves: \t" << boxMoves[3] << std::endl;
            #endif
            for (int k=0; k<4; ++k){
                // std::cout << "K: " << k << std::endl;
                std::bitset<MAXSIZE> curmove = boxMoves[k];
                size_t idx = curmove._Find_first();
                // std::cout << "Direction: " << dirChar[k] << ", Movable boxes:\n";
                while (idx < curmove.size()) {
                    // std::cout << "idx: " << idx << ", (" << idx / totalc << ", " << idx % totalc << ")\n";
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
                    auto [newMovableBoxes, deadFrozen] = getPossibleMoves(newBoxPositions);
                    if (deadFrozen) {
                        idx = curmove._Find_next(idx);
                        continue;
                    }

                    // Compute agent position (where agent must be to push box)
                    int agentRow = boxRow - dr;
                    int agentCol = boxCol - dc;
                    pii newAgentPos = mkp(boxRow, boxCol);
                    pii requiredAgentPos = mkp(agentRow, agentCol);

                    // Update connected component for new agent position
                    auto newCC = connectedComponent(newAgentPos, grid | newBoxPositions);

                    // Update path
                    std::vector<std::pair<pii, int>> newPath = state.path;
                    newPath.push_back({requiredAgentPos, k});

                    State newState = {newAgentPos, newCC, newBoxPositions, newPath};
                    if (isIn(newState, visited)) {
                        idx = curmove._Find_next(idx);
                        continue;
                    }
                    // Check for solution
                    if ((newBoxPositions & ~goalPositions).none()) {
                        std::cout << "States: " << count << std::endl;
                        std::cout << "Moves: " << moveRound << std::endl;
                        std::cout << recoverPath(initialState, newState);
                        std::cout << std::endl;
                        exit(0);
                    }
                    pushin(newState, nv, visited);

                    idx = curmove._Find_next(idx);
                }
            }
            count++;
        }
        v = nv;
        nv.clear();
        moveRound++;
        #ifdef DEBUG
            if (moveRound >= 1) exit(0);
        #endif
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
    std::cout << "totalr: " << totalr << ", totalc: " << totalc << std::endl;
    std::cout << "==========================" << std::endl;
    auto cc = connectedComponent(agentPos, grid | boxPositions);
    // printBitset(goalPositions, {-1, -1});
    simpleDeadlockList();
    // printBitset(simpleDeadStateMask, {-1, -1});
    // std::cout << simpleDeadStateMask << std::endl;
    // exit(0);
    
    bfs({agentPos, cc, boxPositions, {}});

    return 0;
}