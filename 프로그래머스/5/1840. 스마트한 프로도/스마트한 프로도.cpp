#include <vector>
#include <algorithm>

using namespace std;

struct Edge {
    int to;
    int id;
    int type;
};

vector<vector<int>> solution(int n, int m, vector<int> a, vector<int> b, int k, int m1, int m2, vector<int> e1, vector<int> e2) {
    vector<int> state(m + 1, 0);
    for(int x : e1) state[x] += 1;
    for(int x : e2) state[x] += 2;

    vector<vector<Edge>> adj(n + 1);
    for(int i = 1; i <= m; ++i) {
        if(state[i] == 1) {
            adj[a[i-1]].push_back({b[i-1], i, 0});
            adj[b[i-1]].push_back({a[i-1], i, 0});
        } else if(state[i] == 2) {
            adj[a[i-1]].push_back({b[i-1], i, 1});
            adj[b[i-1]].push_back({a[i-1], i, 1});
        }
    }

    vector<bool> visited(m + 1, false);
    vector<vector<pair<int, int>>> paths_11;
    vector<vector<pair<int, int>>> paths_10;
    vector<vector<pair<int, int>>> paths_00;
    vector<vector<pair<int, int>>> cycles;

    auto get_path = [&](int start) {
        vector<pair<int, int>> path;
        int curr = start;
        while(true) {
            bool moved = false;
            for(auto& edge : adj[curr]) {
                if(!visited[edge.id]) {
                    visited[edge.id] = true;
                    path.push_back({edge.id, edge.type});
                    curr = edge.to;
                    moved = true;
                    break;
                }
            }
            if(!moved) break;
        }
        return path;
    };

    for(int i = 1; i <= n; ++i) {
        if(adj[i].size() == 1) {
            bool has_unvisited = false;
            for(auto& e : adj[i]) {
                if(!visited[e.id]) { has_unvisited = true; break; }
            }
            if(has_unvisited) {
                vector<pair<int, int>> p = get_path(i);
                if(!p.empty()) {
                    if(p.front().second == 0 && p.back().second == 1) {
                        reverse(p.begin(), p.end());
                    }
                    if(p.front().second == 1 && p.back().second == 1) {
                        paths_11.push_back(p);
                    } else if(p.front().second == 1 && p.back().second == 0) {
                        paths_10.push_back(p);
                    } else if(p.front().second == 0 && p.back().second == 0) {
                        paths_00.push_back(p);
                    }
                }
            }
        }
    }

    for(int i = 1; i <= n; ++i) {
        if(adj[i].size() == 2) {
            bool has_unvisited = false;
            for(auto& e : adj[i]) {
                if(!visited[e.id]) { has_unvisited = true; break; }
            }
            if(has_unvisited) {
                vector<pair<int, int>> c = get_path(i);
                if(!c.empty()) {
                    if(c[0].second == 1) {
                        int idx = 0;
                        while(idx < c.size() && c[idx].second == 1) idx++;
                        rotate(c.begin(), c.begin() + idx, c.end());
                    }
                    cycles.push_back(c);
                }
            }
        }
    }

    vector<vector<int>> ans;

    auto process_11 = [&](const vector<pair<int, int>>& p) {
        for(size_t i = 0; i + 1 < p.size(); i += 2) {
            ans.push_back({0, p[i+1].first});
            ans.push_back({1, p[i].first});
        }
        ans.push_back({1, p.back().first});
    };

    auto process_10 = [&](const vector<pair<int, int>>& p) {
        for(size_t i = 0; i < p.size(); i += 2) {
            ans.push_back({0, p[i+1].first});
            ans.push_back({1, p[i].first});
        }
    };

    for(const auto& p : paths_11) {
        process_11(p);
    }

    for(const auto& p : paths_10) {
        process_10(p);
    }

    for(const auto& c : cycles) {
        ans.push_back({0, c[0].first});
        vector<pair<int, int>> p;
        for(int i = (int)c.size() - 1; i >= 1; --i) {
            p.push_back(c[i]);
        }
        process_11(p);
    }

    for(const auto& p : paths_00) {
        ans.push_back({0, p[0].first});
        vector<pair<int, int>> sub(p.begin() + 1, p.end());
        process_10(sub);
    }

    return ans;
}