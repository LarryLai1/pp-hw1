#include <bits/stdc++.h>
#include <fstream>
#include <assert.h>
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
    // int agentGroup; // for future use
    std::set<pii> boxPositions; // (row, col)
    int g; // cost from start to current state
    int h; // heuristic: existing '.' count
    int f; // total cost (g + h)
    std::string path; // path taken to reach this state

    // For unordered_set: boxPositions must match, and agentPos must be "reachable" from other's agentPos via agentGoTo.
    bool operator==(const State& other) const {
        return boxPositions == other.boxPositions && agentPos == other.agentPos;
    }
};

// Hash function for State to use in unordered_set
namespace std {
    template <>
    struct hash<State> {
        size_t operator()(const State& s) const {
            size_t res = 0;
            for (const auto& box : s.boxPositions) {
                res ^= (hash<int>()(box.first) << 2) ^ (hash<int>()(box.second) << 3);
            }
            // Also hash agentPos
            res ^= (hash<int>()(s.agentPos.first) << 5) ^ (hash<int>()(s.agentPos.second) << 7);
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

std::vector<std::string> drawGrid(std::vector<std::string> grid, std::set<pii> boxPositions){
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
    std::cout << "Current Grid State:\n";
    for (const auto& row : newGrid) {
        std::cout << row << std::endl;
    }
    return newGrid;
}

std::string agentGoTo(std::pair<int, int> agentPos, std::pair<int, int> targetPos, std::vector<std::string> curgrid) {
    int n = curgrid.size();
    int m = curgrid[0].size();
    std::vector<std::vector<bool>> visited(n, std::vector<bool>(m, false));
    std::queue<std::tuple<int, int, std::string>> q;
    q.push({agentPos.first, agentPos.second, ""});
    visited[agentPos.first][agentPos.second] = true;

    while (!q.empty()) {
        auto [r, c, path] = q.front(); q.pop();
        // std::cout << "Visiting: (" << r << ", " << c << ") with path: " << path << std::endl;
        if (std::make_pair(r, c) == targetPos) {
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
    return ""; // No path found
}

// 0: reachable, 1: simple dead state
std::vector<std::vector<int>> simpleDeadlockList(std::vector<std::string> grid){
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
std::vector<pii> nonFrozenBoxes(std::set<pii> boxPositions, std::vector<std::string> curgrid, 
    std::vector<std::vector<int>> simpleDeadState) {
    std::vector<pii> nonFrozenBoxes;
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
                        (leftCell == 'x' || leftCell == 'X' || rightCell == 'x' || rightCell == 'X') ||
                        ((simpleDeadState[r][c - 1] == 1 || leftCell == '@') 
                        && (simpleDeadState[r][c + 1] == 1 || rightCell == '@'));
        bool vBlocked = (upCell == '#') || (downCell == '#') ||
                        (upCell == 'x' || upCell == 'X' || downCell == 'x' || downCell == 'X') ||
                        ((simpleDeadState[r - 1][c] == 1 || upCell == '@') 
                        && (simpleDeadState[r + 1][c] == 1 || downCell == '@'));
        // std::cout << "Box at (" << r << ", " << c << ") - hBlocked: " << hBlocked << ", vBlocked: " << vBlocked << std::endl;
        if (!(hBlocked && vBlocked)) {
            nonFrozenBoxes.push_back(box);
        }
    }
    return nonFrozenBoxes;
}

// get move directions for a box
// 0: up, 1: down, 2: left, 3: right
std::vector<int> getBoxMoves(std::pair<int, int> boxPos, std::vector<std::string> curgrid, 
    std::vector<std::vector<int>> simpleDeadState, std::pair<int, int> agentPos, 
    std::vector<std::vector<std::string>> &correspondingMoves) {
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
            if (correspondingMoves[ar][ac].empty()){
                correspondingMoves[ar][ac] = agentGoTo(agentPos, {ar, ac}, curgrid);
                // std::cout << correspondingMoves[ar][ac] << std::endl;
            }
            if (agentCell != '#' && agentCell != 'x' && agentCell != 'X' 
                && (!correspondingMoves[ar][ac].empty() || (ar==agentPos.first && ac==agentPos.second))) { // Check if agent can stand to push
                possibleMoves.push_back(i);
            }
        }
    }
    return possibleMoves;
}

// get all box moves
std::vector<std::vector<int>> getAllBoxMoves(std::vector<pii> movableBoxes, std::vector<std::string> curgrid, 
    std::vector<std::vector<int>> simpleDeadState, std::pair<int, int> agentPos, 
    std::vector<std::vector<std::string>> &correspondingMoves) {
    std::vector<std::vector<int>> allBoxMoves;
    for (const auto& box : movableBoxes) {
        auto moves = getBoxMoves(box, curgrid, simpleDeadState, agentPos, correspondingMoves);
        allBoxMoves.push_back(moves);
    }
    return allBoxMoves;
}

bool isDeadState(std::set<pii> boxPositions, std::vector<std::string> curgrid, 
    std::vector<std::vector<int>> simpleDeadState, std::vector<pii> movableBoxes, 
    std::pair<int, int> agentPos, std::vector<std::vector<std::string>> &correspondingMoves){
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

int heuristicFunction(std::set<pii> boxPositions, std::vector<std::string> grid) {
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
    std::unordered_set<State> visited;

    pq.push({agentPos, boxPositions, 0, 0, 0});
    // find simple dead state
    std::vector<std::vector<int>> simpleDeadState = simpleDeadlockList(grid);
    while (!pq.empty()) {
        auto [agentPos, boxPositions, g, h, f, path] = pq.top();
        pq.pop(); visited.insert({agentPos, boxPositions, g, h, f});
        // draw grid
        auto curgrid = drawGrid(grid, boxPositions);
        std::cout << "g: " << g << ", h: " << h << ", f: " << f << ", path: " << path << std::endl;

        // store corresponding moves for agent to reach each position
        std::vector<std::vector<std::string>> correspondingMoves(curgrid.size(), 
                                std::vector<std::string>(curgrid[0].size(), ""));

        // find non-frozen boxes
        auto movableBoxes = nonFrozenBoxes(boxPositions, curgrid, simpleDeadState);
        auto boxMoves = getAllBoxMoves(movableBoxes, curgrid, simpleDeadState, agentPos, correspondingMoves);
        for (int i=0; i<movableBoxes.size(); ++i){
            auto box = movableBoxes[i];
            auto moves = boxMoves[i];
            for (const auto& move : moves) {
                int nr = box.first + dirs[move].first;
                int nc = box.second + dirs[move].second;
                // TODO: maybe check dead state here?
                auto newBoxPositions = boxPositions;
                // update box position
                newBoxPositions.erase(box);
                newBoxPositions.insert(mkp(nr, nc));
                auto heuristic = heuristicFunction(newBoxPositions, grid);
                if (heuristic == 0) {
                    std::cout << "Solution found!\n";
                    std::cout << "Path: " << path + correspondingMoves[box.first][box.second] << 
                        correspondingMoves[box.first][box.second] << "\n";
                    exit(0);
                }
                State newState = {box, newBoxPositions, g+1, heuristic, g+heuristic+1, 
                    path+correspondingMoves[box.first][box.second]+dirChar[move]};
                if (visited.find(newState) != visited.end()){
                    std::cout << "Already visited this state.\n" << std::endl;
                    continue;
                }
                pq.push(newState);
            }
        }
    }

    return 0;
}