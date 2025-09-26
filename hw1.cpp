#include <bits/stdc++.h>
#include <fstream>
#include <assert.h>
#include <boost/functional/hash.hpp>
#include <signal.h>
#include <unistd.h>

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
    std::pair<int, int> agentPos; // (row, col)
    // int agentGroup; // for comparison
    std::set<pii> boxPositions; // (row, col)
    std::string path; // path taken to reach this state

    bool operator==(const State& other) const {
        return boxPositions == other.boxPositions && agentPos == other.agentPos;
    }
};

struct StatePos{
    std::set<pii> agentPositions; // (agent_pos, movable_boxes)
    std::set<pii> boxPositions; // (row, col)
    // Hash function for StatePos using only boxPositions
    bool operator==(const StatePos& other) const {
        return boxPositions == other.boxPositions;
    }
};

namespace std {
    template <>
    struct hash<StatePos> {
        size_t operator()(const StatePos& s) const {
            std::size_t seed = 0;
            for (const auto& box : s.boxPositions) {
                boost::hash_combine(seed, boost::hash<int>()(box.first));
                boost::hash_combine(seed, boost::hash<int>()(box.second));
            }
            return seed;
        }
    };
    
    template <>
    struct hash<std::set<pii>> {
        size_t operator()(const std::set<pii>& s) const {
            std::size_t seed = 0;
            for (const auto& box : s) {
                boost::hash_combine(seed, boost::hash<int>()(box.first));
                boost::hash_combine(seed, boost::hash<int>()(box.second));
            }
            return seed;
        }
    };
}

// Cache for storing movableBoxes under different boxPositions
std::unordered_map<std::set<pii>, std::pair<std::set<pii>, bool>> movableBoxesCache;
int cacheHit = 0;

// 讀取檔案並回傳為vector<string>
std::tuple<std::vector<std::string>, std::set<pii>, pii> readFileToVector(const std::string& filename) {
    std::ifstream file(filename);
    std::vector<std::string> grid;
    
    assert (file.is_open());
    
    std::string line;   
    // 逐行讀取檔案內容
    while (std::getline(file, line)) {
        grid.push_back(line);
    }
    file.close();
    std::set<pii> boxPositions;
    pii agentPos;
    for (int i = 0; i < grid.size(); ++i) {
        for (int j = 0; j < grid[i].size(); ++j) {
            char cell = grid[i][j];
            switch (cell){
                case 'o':
                    agentPos = mkp(i, j);
                    grid[i][j] = ' ';
                    continue;
                case 'O':
                    agentPos = mkp(i, j);
                    grid[i][j] = '.';
                    continue;
                case '!':
                    agentPos = mkp(i, j);
                    grid[i][j] = '@';
                    continue;
                case 'x':
                    boxPositions.insert(mkp(i, j));
                    grid[i][j] = ' ';
                    continue;
                case 'X':
                    boxPositions.insert(mkp(i, j));
                    grid[i][j] = '.';
                    continue;
                default:
                    continue;
            }
        }
    }
    return std::make_tuple(grid, boxPositions, agentPos);
}

std::vector<std::string> drawGrid(const std::vector<std::string>& grid, const std::set<pii>& boxPositions){
    std::vector<std::string> newGrid = grid;
    for (const auto& box : boxPositions) {
        int r = box.first;
        int c = box.second;
        if (newGrid[r][c] == '.') {
            newGrid[r][c] = 'X'; // Box on target
        } else {
            newGrid[r][c] = 'x'; // Box on regular tile
        }
    }
    return newGrid;
}

std::string agentGoTo(const std::pair<int, int>& agentPos, const std::pair<int, int>& targetPos, 
    const std::vector<std::string>& curgrid) {
    // std::cout << "Finding path for agent from (" << agentPos.first << ", " << agentPos.second << ") to (" << targetPos.first << ", " << targetPos.second << ")\n";
    if (agentPos == targetPos) return "";
    int n = curgrid.size();
    int m = curgrid[0].size();
    int tr = targetPos.first, tc = targetPos.second;
    std::vector<std::vector<bool>> visited(n, std::vector<bool>(m, false));
    std::queue<std::tuple<int, int, std::string>> q;
    q.push({agentPos.first, agentPos.second, ""});
    visited[agentPos.first][agentPos.second] = true;
    while (!q.empty()) {
        auto [r, c, path] = q.front(); q.pop();
        // std::cout << "Visiting: (" << r << ", " << c << ") from (" << agentPos.first << ", " << agentPos.second << ") with path: " << path << std::endl;
        if (r == tr && c == tc) {
            // std::cout << "Path found: " << path << std::endl;
            return path;
        }
        for (int d = 0; d < 4; ++d) {
            int nr = r + dirs[d].first;
            int nc = c + dirs[d].second;
            if (nr < 0 || nr >= n || nc < 0 || nc >= m) continue;
            char cell = curgrid[nr][nc];
            if (cell == '#' || cell == 'x' || cell == 'X') continue;
            if (!visited[nr][nc]) {
                visited[nr][nc] = true;
                q.push({nr, nc, path + dirChar[d]});
            }
        }
    }
    return "N"; // No path found
}

bool hasPath(const pii& agentPos, const pii& targetPos, const std::set<pii>& boxPositions, 
    const std::vector<std::string>& grid) {
    int n = grid.size();
    int m = grid[0].size();
    int tr = targetPos.first, tc = targetPos.second;
    std::vector<std::vector<bool>> visited(n, std::vector<bool>(m, false));
    std::queue<pii> q;
    q.push(agentPos);
    visited[agentPos.first][agentPos.second] = true;
    while (!q.empty()) {
        auto [r, c] = q.front(); q.pop();
        if (r == tr && c == tc) {
            return true;
        }
        for (int d = 0; d < 4; ++d) {
            int nr = r + dirs[d].first;
            int nc = c + dirs[d].second;
            if (nr < 0 || nr >= n || nc < 0 || nc >= m) continue;
            char cell = grid[nr][nc];
            if (cell == '#' || boxPositions.count({nr, nc})) continue;
            if (!visited[nr][nc]) {
                visited[nr][nc] = true;
                q.push({nr, nc});
            }
        }
    }
    return false; // No path found
}

// 0: reachable, 1: simple dead state
std::vector<std::vector<int>> simpleDeadlockList(const std::vector<std::string>& grid){
    std::queue<std::pair<int, int>> goalPositions;
    std::vector<std::vector<int>> simpleDeadState(grid.size(), std::vector<int>(grid[0].size(), 0));
    for (int i = 0; i < grid.size(); ++i) {
        for (int j = 0; j < grid[i].size(); ++j) {
            char cell = grid[i][j];
            if (cell == '#') continue; // wall
            if (cell == '.') {
                goalPositions.push(mkp(i, j));
                continue; // target
            }
            // empty space
            simpleDeadState[i][j] = 1;
        }
    }

    while (goalPositions.size() > 0) {
        std::pair<int, int> goal = goalPositions.front(); goalPositions.pop();
        int r = goal.first;
        int c = goal.second;
        // Check all four directions
        for (const auto& dir : dirs) {
            int nr = r + dir.first;
            int nc = c + dir.second;
            if (simpleDeadState[nr][nc] == 1 && grid[nr+dir.first][nc+dir.second] != '#') {
                simpleDeadState[nr][nc] = 0; // Mark as reachable
                goalPositions.push(mkp(nr, nc));
            }
        }
    }
    return simpleDeadState;
}


// return all "currently" non-frozen boxes
std::pair<std::set<pii>, bool> getMovableBoxes(const std::set<pii>& boxPositions, const std::vector<std::string>& curgrid, 
    const std::vector<std::vector<int>>& simpleDeadState) {
    
    // Check cache first
    auto cacheIt = movableBoxesCache.find(boxPositions);
    if (cacheIt != movableBoxesCache.end()) {
        cacheHit++;
        return movableBoxesCache[boxPositions];
    }
    
    std::set<pii> movableBoxes;
    std::vector<pii> nonFrozenBoxes, FrozenBoxes, tempBoxes;
    std::vector<std::vector<int>> frozen(curgrid.size(), std::vector<int>(curgrid[0].size(), -1));
    for (const auto& box : boxPositions) {
        int r = box.first;
        int c = box.second;
        // std::cout << "Analyzing box at: (" << r << ", " << c << ")\n";
        
        char leftCell = curgrid[r][c - 1];
        char rightCell = curgrid[r][c + 1];
        char upCell = curgrid[r - 1][c];
        char downCell = curgrid[r + 1][c];

        // Check if frozen under no box condition
        bool hBlocked = (leftCell == '#' || rightCell == '#' ) ||
                        ((simpleDeadState[r][c - 1] == 1 || leftCell == '@') 
                        && (simpleDeadState[r][c + 1] == 1 || rightCell == '@'));
        bool vBlocked = (upCell == '#') || (downCell == '#') ||
                        ((simpleDeadState[r - 1][c] == 1 || upCell == '@') 
                        && (simpleDeadState[r + 1][c] == 1 || downCell == '@'));
        bool hBoxed = (leftCell == 'x' || leftCell == 'X' || rightCell == 'x' || rightCell == 'X');
        bool vBoxed = (upCell == 'x' || upCell == 'X' || downCell == 'x' || downCell == 'X');
        // std::cout << "Box at (" << r << ", " << c << ") - hBlocked: " << hBlocked << ", vBlocked: " << vBlocked << std::endl;
        if (!((hBlocked || hBoxed) && (vBlocked || vBoxed))) {
            movableBoxes.insert(box);
            frozen[r][c] = 0;
        }
        else if (hBlocked && vBlocked){
            frozen[r][c] = 1;
        } else {
            FrozenBoxes.push_back(box);
            frozen[r][c] = 1;
        }
    }
    bool changed = true;
    while (!FrozenBoxes.empty() && changed){
        changed = false;
        for (const auto &cell: FrozenBoxes) {
            auto [r, c] = cell;
            char leftCell = curgrid[r][c - 1];
            char rightCell = curgrid[r][c + 1];
            char upCell = curgrid[r - 1][c];
            char downCell = curgrid[r + 1][c];

            bool hBlocked = (leftCell == '#' || rightCell == '#' ) ||
                            (frozen[r][c-1] > 0 || frozen[r][c+1] > 0) ||
                            ((simpleDeadState[r][c - 1] == 1 || leftCell == '@') 
                            && (simpleDeadState[r][c + 1] == 1 || rightCell == '@'));
            bool vBlocked = (upCell == '#') || (downCell == '#') ||
                            (frozen[r-1][c] > 0 || frozen[r+1][c] > 0) ||
                            ((simpleDeadState[r - 1][c] == 1 || upCell == '@') 
                            && (simpleDeadState[r + 1][c] == 1 || downCell == '@'));
            if (!(hBlocked && vBlocked)) {
                nonFrozenBoxes.push_back(cell);
                changed = true;
                frozen[r][c] = 0;
            } else {
                tempBoxes.push_back(cell);
            }
        }
        FrozenBoxes = tempBoxes;
        tempBoxes.clear();
    }
    bool deadFrozen=false;
    for (const auto &cell: FrozenBoxes) {
        auto [r, c] = cell;
        if (curgrid[r][c] != 'X') {
            deadFrozen = true;
        }
    }
    // Store result in cache before returning
    std::pair<std::set<pii>, bool> result = {movableBoxes, deadFrozen};
    movableBoxesCache[boxPositions] = result;
    return result;
}

// get move directions for a box
// 0: up, 1: down, 2: left, 3: right
std::vector<int> getBoxMoves(const std::pair<int, int>& boxPos, const std::vector<std::string>& curgrid, 
    const std::vector<std::vector<int>>& simpleDeadState, const std::pair<int, int>& agentPos, 
    std::vector<std::vector<std::string>>& correspondingMoves) {
    std::vector<int> possibleMoves;
    int r = boxPos.first;
    int c = boxPos.second;
    // Check all four directions for possible pushes
    for (int i = 0; i < 4; ++i) {
        auto dir = dirs[i];
        int nr = r + dir.first; // new row for box
        int nc = c + dir.second; // new col for box
        int ar = r - dir.first; // agent row to push
        int ac = c - dir.second; // agent col to push
        char destCell = curgrid[nr][nc], agentCell = curgrid[ar][ac];
        if ((destCell == ' ' || destCell == '.') && (simpleDeadState[nr][nc] == 0)) { // Check if box can be moved to new position
            if (correspondingMoves[ar][ac] == "N"){ // Not calculated yet
                correspondingMoves[ar][ac] = agentGoTo(agentPos, {ar, ac}, curgrid);
                // std::cout << correspondingMoves[ar][ac] << std::endl;
            }
            if (agentCell != '#' && agentCell != 'x' && agentCell != 'X' 
                && (correspondingMoves[ar][ac] != "N" || (ar==agentPos.first && ac==agentPos.second))) { // Check if agent can stand to push
                possibleMoves.push_back(i);
            }
        }
    }
    return possibleMoves;
}

// get all box moves
std::vector<std::pair<pii, std::vector<int>>> getAllBoxMoves(const std::set<pii>& movableBoxes, const std::vector<std::string>& curgrid, 
    const std::vector<std::vector<int>>& simpleDeadState, const std::pair<int, int>& agentPos, 
    std::vector<std::vector<std::string>>& correspondingMoves) {
    std::vector<std::pair<pii, std::vector<int>>> allBoxMoves;
    for (const auto& box : movableBoxes) {
        auto moves = getBoxMoves(box, curgrid, simpleDeadState, agentPos, correspondingMoves);
        allBoxMoves.push_back({box, moves});
    }
    return allBoxMoves;
}

int oldHeuristicFunction(const std::set<pii>& boxPositions, const std::vector<std::string>& grid) {
    int h = 0;
    for (const auto& box : boxPositions) {
        auto [r, c] = box;
        if (grid[r][c] != '.') h++;
    }
    return h;
}

bool isIn(const State& state, const std::unordered_set<StatePos> &visited,
    const std::vector<std::string>& grid){
    StatePos statePos;
    statePos.agentPositions.insert(state.agentPos);
    statePos.boxPositions = state.boxPositions;
    auto found = visited.find(statePos);
    if (found != visited.end()) {
        for (const auto& pos : found->agentPositions) {
            // TODO: duplication check can be improved
            // reachable from previous position
            if (hasPath(pos, state.agentPos, state.boxPositions, grid)) {
                return true;
            }
        }
    }
    return false;
}

void pushin(const State& state, std::vector<State> &v, std::unordered_set<StatePos> &visited,
    const std::vector<std::string>& grid) {
    StatePos statePos;
    statePos.agentPositions.insert(state.agentPos);
    statePos.boxPositions = state.boxPositions;
    auto found = visited.find(statePos);
    if (found != visited.end()) {
        for (const auto& pos : found->agentPositions) {
            // TODO: duplication check can be improved
            // reachable from previous position
            if (hasPath(pos, state.agentPos, state.boxPositions, grid)) {
                return;
            }
        }
        // not reachable, add new agent position
        statePos.agentPositions = found->agentPositions;
        // std::cout << "A\n";
        statePos.agentPositions.insert(state.agentPos);
        // std::cout << "B\n";
        visited.insert(statePos);
        // std::cout << "C\n";
        v.push_back(state);
        // std::cout << "D\n";
        visited.erase(found);
        // std::cout << "E\n";
        return;
    }
    // not found, insert new state
    visited.insert(statePos);
    v.push_back(state);
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
    auto [grid, boxPositions, agentPos] = readFileToVector(filename);
    std::vector<State> v, nv;
    std::unordered_set<StatePos> visited;

    // Initialize with movable boxes for initial state
    auto curgrid = drawGrid(grid, boxPositions);
    std::vector<std::vector<int>> simpleDeadState = simpleDeadlockList(grid);
    auto [initialMovableBoxes, initialDeadFrozen] = getMovableBoxes(boxPositions, curgrid, simpleDeadState);
    if (!initialDeadFrozen) {
        pushin({agentPos, boxPositions, ""}, v, visited, grid);
    }
    
    int count = 0;
    int moveRound = 0;
    while (!v.empty()) {
        for (const auto& state: v){
            auto [agentPos, boxPositions, path] = state;
            auto movableBoxes = movableBoxesCache[boxPositions].first;
            // draw grid
            auto curgrid = drawGrid(grid, boxPositions);
    
            // store corresponding moves for agent to reach each position
            std::vector<std::vector<std::string>> correspondingMoves(curgrid.size(), 
                                    std::vector<std::string>(curgrid[0].size(), "N"));
    
            auto boxMoves = getAllBoxMoves(movableBoxes, curgrid, simpleDeadState, agentPos, correspondingMoves);
            for (const auto& [box, moves] : boxMoves){
                for (const auto& move : moves) {
                    int nr = box.first + dirs[move].first;
                    int nc = box.second + dirs[move].second;
                    // TODO: maybe check dead state here?
                    auto newBoxPositions = boxPositions;
                    // update box position
                    newBoxPositions.erase(box);
                    newBoxPositions.insert(mkp(nr, nc));
                    
                    // Draw new grid and get movable boxes before pushing
                    auto newCurgrid = drawGrid(grid, newBoxPositions);
                    auto [newMovableBoxes, deadFrozen] = getMovableBoxes(newBoxPositions, newCurgrid, simpleDeadState);
                    if (deadFrozen) {
                        continue;
                    }
                    
                    auto heuristic = oldHeuristicFunction(newBoxPositions, grid);
                    State newState = {box, newBoxPositions, 
                        path+correspondingMoves[box.first - dirs[move].first][box.second - dirs[move].second]+dirChar[move]};
                    if (isIn(newState, visited, grid)) {
                        continue;
                    }
                    if (heuristic == 0) {
                        std::cout << "States: " << count << std::endl;
                        std::cout << "Moves: " << moveRound << std::endl;
                        std::cout << "Cache hit: " << cacheHit << std::endl;
                        std::cout << newState.path << std::endl;
                        exit(0);
                    }
                    pushin(newState, nv, visited, grid);
                }
            }
            count++;
        }
        v = nv;
        nv.clear();
        moveRound++;
    }
    return 0;
}