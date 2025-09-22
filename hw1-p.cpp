#include <bits/stdc++.h>
#include <fstream>
#include <assert.h>
#include <omp.h>
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
    int g; // cost from start to current state
    int h; // heuristic: existing '.' count
    int f; // total cost (g + h)
    std::string path; // path taken to reach this state

    bool operator==(const State& other) const {
        return boxPositions == other.boxPositions && agentPos == other.agentPos;
    }
};

struct StatePos{
    std::set<pii> agentPositions; // (row, col)
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
            size_t res = 0;
            for (const auto& box : s.boxPositions) {
                res ^= (hash<int>()(box.first) << 2) ^ (hash<int>()(box.second) << 3);
            }
            return res;
        }
    };
}


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
    #ifdef DEBUG
    std::cout << "Current Grid State:\n";
    for (const auto& row : newGrid) {
        std::cout << row << std::endl;
    }
    #endif
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
    // Print simpleDeadState for debugging
    #ifdef DEBUG
    std::cout << "Simple Dead State Map:\n";
    for (int i = 0; i < simpleDeadState.size(); ++i) {
        for (int j = 0; j < simpleDeadState[i].size(); ++j) {
            if (simpleDeadState[i][j] == 1)
                std::cout << "1";
            else
            std::cout << grid[i][j];
        }
        std::cout << std::endl;
    }
    #endif
    return simpleDeadState;
}

// return all "currently" non-frozen boxes
std::pair<std::vector<pii>, bool> getMovableBoxes(const std::set<pii>& boxPositions, const std::vector<std::string>& curgrid, 
    const std::vector<std::vector<int>>& simpleDeadState) {
    std::vector<pii> movableBoxes, nonFrozenBoxes, FrozenBoxes, tempBoxes;
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
            movableBoxes.push_back(box);
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
            #ifdef DEBUG
            std::cout << "Box at (" << r << ", " << c << ") is dead frozen.\n";
            #endif
        }
    }
    return {movableBoxes, deadFrozen};
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
std::vector<std::vector<int>> getAllBoxMoves(const std::vector<pii>& movableBoxes, const std::vector<std::string>& curgrid, 
    const std::vector<std::vector<int>>& simpleDeadState, const std::pair<int, int>& agentPos, 
    std::vector<std::vector<std::string>>& correspondingMoves) {
    std::vector<std::vector<int>> allBoxMoves;
    for (const auto& box : movableBoxes) {
        auto moves = getBoxMoves(box, curgrid, simpleDeadState, agentPos, correspondingMoves);
        allBoxMoves.push_back(moves);
    }
    return allBoxMoves;
}

bool isDeadState(const std::set<pii>& boxPositions, const std::vector<std::string>& curgrid, 
    const std::vector<std::vector<int>>& simpleDeadState, const std::vector<pii>& movableBoxes, 
    const std::pair<int, int>& agentPos, std::vector<std::vector<std::string>>& correspondingMoves) {
    /*
    1. check if any box is in a simple dead state
    2. check if there is no movable box
    3. check if all movable boxes have no possible moves
    */
    for (const auto& box : boxPositions) {
        int r = box.first;
        int c = box.second;
        if (simpleDeadState[r][c] == 1) {
            // std::cout << "Box at (" << r << ", " << c << ") is in a simple dead state.\n";
            return true;
        }
    }
    if (movableBoxes.empty()) {
        // std::cout << "No movable boxes available.\n";
        return true;
    }
    for (const auto& box : movableBoxes) {
        auto possibleMoves = getBoxMoves(box, curgrid, simpleDeadState, agentPos, correspondingMoves);
        if (!possibleMoves.empty()) {
            return false; // Found a movable box with possible moves
        }
    }
    return true; // All checks passed, it's a dead state
}

int heuristicFunction(const std::set<pii>& boxPositions, const std::vector<std::string>& grid) {
    int h = 0;
    for (const auto& box : boxPositions) {
        int r = box.first;
        int c = box.second;
        if (grid[r][c] != '.') {
            h += 1; // Box on regular tile
        }
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
        #ifdef DEBUG
        std::cout << "Found existing state, checking reachability...\n";
        std::cout << "Box positions: ";
        for (const auto& box : state.boxPositions) {
            std::cout << "(" << box.first << ", " << box.second << ") ";
        }
        std::cout << "\n";
        #endif
        for (const auto& pos : found->agentPositions) {
            // TODO: duplication check can be improved
            // reachable from previous position
            #ifdef DEBUG
            std::cout << "Checking if agent can reach from (" << pos.first << ", " << pos.second << ") to (" 
                      << state.agentPos.first << ", " << state.agentPos.second << ")\n";
            #endif
            if (hasPath(pos, state.agentPos, state.boxPositions, grid)) {
                #ifdef DEBUG
                std::cout << "Agent can reach the new position.\n";
                #endif
                return true;
            }
            else {
                #ifdef DEBUG
                std::cout << "Agent cannot reach the new position from (" << pos.first << ", " << pos.second << ").\n";
                #endif
            }
        }
    }
    #ifdef DEBUG
    std::cout << "No existing state found or no reachable agent position.\n";
    #endif
    return false;
}

void pushin(const State& state, std::priority_queue<State, std::vector<State>, 
    std::function<bool(const State&, const State&)>> &pq, std::unordered_set<StatePos> &visited,
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
                #ifdef DEBUG
                std::cout << "Agent can reach the new position, not pushing new state.\n";
                #endif
                return;
            }
        }
        // not reachable, add new agent position
        #ifdef DEBUG
        std::cout << "Agent cannot reach the new position, adding new agent position to existing state.\n";
        #endif
        statePos.agentPositions = found->agentPositions;
        // std::cout << "A\n";
        statePos.agentPositions.insert(state.agentPos);
        // std::cout << "B\n";
        visited.insert(statePos);
        // std::cout << "C\n";
        pq.push(state);
        // std::cout << "D\n";
        visited.erase(found);
        // std::cout << "E\n";
        #ifdef DEBUG
        std::cout << "Done\n";
        #endif
        return;
    }
    // not found, insert new state
    #ifdef DEBUG
    std::cout << "Inserting new state into visited set and priority queue.\n";
    #endif
    visited.insert(statePos);
    pq.push(state);
    #ifdef DEBUG
    std::cout << "Done\n";
    #endif
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
    std::priority_queue<State, std::vector<State>, 
        std::function<bool(const State&, const State&)>> pq(
        [](const State& a, const State& b) { return a.f > b.f; });
    std::unordered_set<StatePos> visited;

    pushin({agentPos, boxPositions, 0, 0, 0, ""}, pq, visited, grid);
    // find simple dead state
    std::vector<std::vector<int>> simpleDeadState = simpleDeadlockList(grid);
    while (!pq.empty()) {
        auto [agentPos, boxPositions, g, h, f, path] = pq.top();
        #ifdef DEBUG
        std::cout << "-------------------------------------------\n";
        std::cout << "Exploring state with g: " << g << ", h: " << h << ", f: " << f << ", path: " << path << std::endl;
        std::cout << "Agent Position: (" << agentPos.first << ", " << agentPos.second << ")\n";
        #endif
        pq.pop();
        // draw grid
        auto curgrid = drawGrid(grid, boxPositions);
        // std::cout << "g: " << g << ", h: " << h << ", f: " << f << ", path: " << path << std::endl;

        // store corresponding moves for agent to reach each position
        std::vector<std::vector<std::string>> correspondingMoves(curgrid.size(), 
                                std::vector<std::string>(curgrid[0].size(), "N"));

        // find non-frozen boxes
        auto [movableBoxes, deadFrozen] = getMovableBoxes(boxPositions, curgrid, simpleDeadState);
        if (deadFrozen) {
            #ifdef DEBUG
            std::cout << "Dead frozen box detected. Skipping this state.\n";
            #endif
            continue;
        }
        auto boxMoves = getAllBoxMoves(movableBoxes, curgrid, simpleDeadState, agentPos, correspondingMoves);
        for (int i=0; i<movableBoxes.size(); ++i){
            auto box = movableBoxes[i];
            auto moves = boxMoves[i];
            for (const auto& move : moves) {
                int nr = box.first + dirs[move].first;
                int nc = box.second + dirs[move].second;
                #ifdef DEBUG
                std::cout << "** Checking Move " << dirChar[move] << " for Box at (" << box.first << ", " << box.second << ") to (" << nr << ", " << nc << ")\n";
                #endif
                // TODO: maybe check dead state here?
                auto newBoxPositions = boxPositions;
                // update box position
                newBoxPositions.erase(box);
                newBoxPositions.insert(mkp(nr, nc));
                auto heuristic = heuristicFunction(newBoxPositions, grid);
                State newState = {box, newBoxPositions, g+1, heuristic, g+heuristic+1, 
                    path+correspondingMoves[box.first - dirs[move].first][box.second - dirs[move].second]+dirChar[move]};
                if (isIn(newState, visited, grid)) {
                    #ifdef DEBUG
                    std::cout << "Already existing this state." << std::endl;
                    #endif
                    continue;
                }
                if (heuristic == 0) {
                    #ifdef DEBUG
                    std::cout << "Solution found!\n";
                    std::cout << "Length: " << newState.g << std::endl;
                    std::cout << "Path: ";
                    #endif
                    std::cout << newState.path << std::endl;
                    exit(0);
                }
                #ifdef DEBUG
                std::cout << "Box: (" << box.first << ", " << box.second << "), Move: " << dirChar[move] << ", Push path: " << newState.path << std::endl;
                #endif
                pushin(newState, pq, visited, grid);
            }
        }
    }

    return 0;
}