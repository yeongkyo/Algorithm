#include <iostream>
#include <vector>
#include <string>
#include <queue>
#include <map>

using namespace std;

struct Point {
    int r, c;
    bool operator<(const Point& o) const {
        if (r != o.r) return r < o.r;
        return c < o.c;
    }
    bool operator==(const Point& o) const {
        return r == o.r && c == o.c;
    }
};

struct Edge {
    Point p1, p2;
    Edge(Point a, Point b) {
        if (b < a) { p1 = b; p2 = a; }
        else { p1 = a; p2 = b; }
    }
    bool operator<(const Edge& o) const {
        if (!(p1 == o.p1)) return p1 < o.p1;
        return p2 < o.p2;
    }
};

int get_type(int r, int c, const vector<string>& grid) {
    int hash_count = 0;
    for(int i=0; i<3; ++i) {
        for(int j=0; j<3; ++j) {
            if(grid[3*r+i][3*c+j] == '#') hash_count++;
        }
    }
    if (hash_count == 0) return 0;
    if (hash_count == 9) return 1;
    if (hash_count == 8) {
        char center = grid[3*r+1][3*c+1];
        if (center >= '0' && center <= '4') return 10 + (center - '0');
    }
    if (hash_count == 6) {
        if (grid[3*r+2][3*c+2] == '.') return 2;
        if (grid[3*r+2][3*c+0] == '.') return 3;
        if (grid[3*r+0][3*c+2] == '.') return 4;
        if (grid[3*r+0][3*c+0] == '.') return 5;
    }
    return -1;
}

vector<Edge> get_edges(int type, int r, int c) {
    vector<Edge> edges;
    Point TL = {r, c};
    Point TR = {r, c+1};
    Point BL = {r+1, c};
    Point BR = {r+1, c+1};
    if (type == 0) {
        edges.push_back(Edge(TL, TR));
        edges.push_back(Edge(TR, BR));
        edges.push_back(Edge(BR, BL));
        edges.push_back(Edge(BL, TL));
    } else if (type == 2) {
        edges.push_back(Edge(TR, BR));
        edges.push_back(Edge(BR, BL));
        edges.push_back(Edge(BL, TR));
    } else if (type == 3) {
        edges.push_back(Edge(BL, TL));
        edges.push_back(Edge(BR, BL));
        edges.push_back(Edge(TL, BR));
    } else if (type == 4) {
        edges.push_back(Edge(TL, TR));
        edges.push_back(Edge(TR, BR));
        edges.push_back(Edge(TL, BR));
    } else if (type == 5) {
        edges.push_back(Edge(TL, TR));
        edges.push_back(Edge(BL, TL));
        edges.push_back(Edge(BL, TR));
    }
    return edges;
}

bool has_right_white(int type) { return type == 0 || type == 2 || type == 4; }
bool has_left_white(int type) { return type == 0 || type == 3 || type == 5; }
bool has_bottom_white(int type) { return type == 0 || type == 2 || type == 3; }
bool has_top_white(int type) { return type == 0 || type == 4 || type == 5; }

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    int N, M;
    if (!(cin >> N >> M)) return 0;
    vector<string> grid(3 * N);
    for (int i = 0; i < 3 * N; ++i) {
        cin >> grid[i];
    }
    
    vector<vector<int>> blocks(N, vector<int>(M));
    for (int r = 0; r < N; ++r) {
        for (int c = 0; c < M; ++c) {
            blocks[r][c] = get_type(r, c, grid);
        }
    }

    int dr[] = {-1, 1, 0, 0};
    int dc[] = {0, 0, -1, 1};
    for (int r = 0; r < N; ++r) {
        for (int c = 0; c < M; ++c) {
            if (blocks[r][c] >= 10) {
                int expected = blocks[r][c] - 10;
                int count = 0;
                for (int i = 0; i < 4; ++i) {
                    int nr = r + dr[i];
                    int nc = c + dc[i];
                    if (nr >= 0 && nr < N && nc >= 0 && nc < M) {
                        int t = blocks[nr][nc];
                        if (t >= 2 && t <= 5) count++;
                    }
                }
                if (count != expected) {
                    cout << "NO\n";
                    return 0;
                }
            }
        }
    }

    vector<vector<int>> adj_blocks(N * M);
    for (int r = 0; r < N; ++r) {
        for (int c = 0; c < M; ++c) {
            int u = r * M + c;
            int type_u = blocks[r][c];
            if (type_u == 1 || type_u >= 10) continue;

            if (c + 1 < M) {
                int type_v = blocks[r][c+1];
                if (type_v != 1 && type_v < 10) {
                    if (has_right_white(type_u) && has_left_white(type_v)) {
                        adj_blocks[u].push_back(r * M + c + 1);
                        adj_blocks[r * M + c + 1].push_back(u);
                    }
                }
            }
            if (r + 1 < N) {
                int type_v = blocks[r+1][c];
                if (type_v != 1 && type_v < 10) {
                    if (has_bottom_white(type_u) && has_top_white(type_v)) {
                        adj_blocks[u].push_back((r + 1) * M + c);
                        adj_blocks[(r + 1) * M + c].push_back(u);
                    }
                }
            }
        }
    }

    vector<bool> visited(N * M, false);
    for (int i = 0; i < N * M; ++i) {
        int r = i / M;
        int c = i % M;
        if (blocks[r][c] == 1 || blocks[r][c] >= 10) continue;
        if (visited[i]) continue;

        vector<int> comp;
        queue<int> q;
        q.push(i);
        visited[i] = true;
        while(!q.empty()) {
            int u = q.front(); q.pop();
            comp.push_back(u);
            for (int v : adj_blocks[u]) {
                if (!visited[v]) {
                    visited[v] = true;
                    q.push(v);
                }
            }
        }

        map<Edge, int> edge_counts;
        for (int u : comp) {
            int ur = u / M;
            int uc = u % M;
            vector<Edge> edges = get_edges(blocks[ur][uc], ur, uc);
            for (Edge e : edges) {
                edge_counts[e]++;
            }
        }

        vector<Edge> boundary;
        for (auto const& kv : edge_counts) {
            if (kv.second == 1) boundary.push_back(kv.first);
            else if (kv.second > 2) {
                cout << "NO\n";
                return 0;
            }
        }

        if (boundary.empty()) continue;

        map<Point, vector<Point>> poly_adj;
        for (Edge e : boundary) {
            poly_adj[e.p1].push_back(e.p2);
            poly_adj[e.p2].push_back(e.p1);
        }

        for (auto const& kv : poly_adj) {
            if (kv.second.size() != 2) {
                cout << "NO\n";
                return 0;
            }
        }

        Point start_pt = poly_adj.begin()->first;
        vector<Point> loop;
        Point curr = start_pt;
        Point prev = {-1, -1};
        do {
            loop.push_back(curr);
            Point next_pt = poly_adj[curr][0];
            if (next_pt == prev) next_pt = poly_adj[curr][1];
            prev = curr;
            curr = next_pt;
        } while (!(curr == start_pt));

        if (loop.size() != boundary.size()) {
            cout << "NO\n";
            return 0;
        }

        int left_turns = 0, right_turns = 0;
        int L = loop.size();
        for (int j = 0; j < L; ++j) {
            Point p1 = loop[(j - 1 + L) % L];
            Point p2 = loop[j];
            Point p3 = loop[(j + 1) % L];
            
            long long dr1 = p2.r - p1.r;
            long long dc1 = p2.c - p1.c;
            long long dr2 = p3.r - p2.r;
            long long dc2 = p3.c - p2.c;

            long long Z = dr1 * dc2 - dc1 * dr2;
            long long D = dr1 * dr2 + dc1 * dc2;

            if (Z == 0 && D > 0) {
                continue;
            } else if (D == 0 && Z != 0) {
                if (Z > 0) left_turns++;
                else right_turns++;
            } else {
                cout << "NO\n";
                return 0;
            }
        }

        if (!((left_turns == 4 && right_turns == 0) || (left_turns == 0 && right_turns == 4))) {
            cout << "NO\n";
            return 0;
        }
    }

    cout << "YES\n";
    return 0;
}