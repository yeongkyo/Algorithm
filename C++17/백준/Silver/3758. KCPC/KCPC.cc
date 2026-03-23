#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

void solve() {
    int n, k, t, m;
    cin >> n >> k >> t >> m;

    vector<vector<int>> scores(n + 1, vector<int>(k + 1, 0));
    vector<int> submit_cnt(n + 1, 0);
    vector<int> last_submit(n + 1, 0);

    for (int time_idx = 1; time_idx <= m; time_idx++) {
        int id, p, s;
        cin >> id >> p >> s;
        
        scores[id][p] = max(scores[id][p], s);
        submit_cnt[id]++;
        last_submit[id] = time_idx;
    }

    vector<int> total_scores(n + 1, 0);
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= k; j++) {
            total_scores[i] += scores[i][j];
        }
    }

    int my_score = total_scores[t];
    int my_cnt = submit_cnt[t];
    int my_time = last_submit[t];

    int rank = 1;
    for (int i = 1; i <= n; i++) {
        if (i == t) continue;

        if (total_scores[i] > my_score) {
            rank++;
        } else if (total_scores[i] == my_score) {
            if (submit_cnt[i] < my_cnt) {
                rank++;
            } else if (submit_cnt[i] == my_cnt) {
                if (last_submit[i] < my_time) {
                    rank++;
                }
            }
        }
    }
    cout << rank << "\n";
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    int T;
    if (cin >> T) {
        while (T--) {
            solve();
        }
    }
    return 0;
}