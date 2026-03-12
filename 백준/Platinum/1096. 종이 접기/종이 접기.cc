#include <iostream>
#include <vector>
#include <algorithm>
#include <set>
#include <queue>

using namespace std;
typedef vector<int> Pile;
typedef vector<Pile> State;

State normalize(State s) {
    for (auto& p : s) sort(p.begin(), p.end());
    State rev = s;
    reverse(rev.begin(), rev.end());
    return (rev < s) ? rev : s;
}

vector<Pile> get_all_piles(int L) {
    State initial;
    for (int i = 0; i < L; ++i) initial.push_back({i});

    set<State> seen;
    set<Pile> unique_piles;
    queue<State> q;

    State start = normalize(initial);
    seen.insert(start);
    q.push(start);

    while (!q.empty()) {
        State curr = q.front();
        q.pop();

        for (const auto& p : curr) unique_piles.insert(p);

        int n = curr.size();
        for (int k = 1; k < n; ++k) {
            int next_n = max(k, n - k);
            State next_s(next_n);
            for (int i = 0; i < next_n; ++i) {
                int left = k - 1 - i;
                int right = k + i;
                if (left >= 0) {
                    for (int v : curr[left]) next_s[i].push_back(v);
                }
                if (right < n) {
                    for (int v : curr[right]) next_s[i].push_back(v);
                }
            }
            State norm_next = normalize(next_s);
            if (seen.find(norm_next) == seen.end()) {
                seen.insert(norm_next);
                q.push(norm_next);
            }
        }
    }
    return vector<Pile>(unique_piles.begin(), unique_piles.end());
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int N, M;
    if (!(cin >> N >> M)) return 0;

    vector<vector<int>> grid(N, vector<int>(M));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            cin >> grid[i][j];
        }
    }

    vector<Pile> row_sets = get_all_piles(N);
    vector<Pile> col_sets = get_all_piles(M);

    long long max_sum = -2e18;

    for (const auto& rs : row_sets) {
        for (const auto& cs : col_sets) {
            long long current_sum = 0;
            for (int r : rs) {
                for (int c : cs) {
                    current_sum += grid[r][c];
                }
            }
            if (current_sum > max_sum) max_sum = current_sum;
        }
    }

    cout << max_sum << endl;

    return 0;
}