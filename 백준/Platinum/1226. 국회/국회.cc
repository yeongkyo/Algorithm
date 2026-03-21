#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>

using namespace std;

struct Party {
    int seats;
    int id;
};

bool dp[305][50005];

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    if (!(cin >> n)) return 0;

    vector<Party> parties(n + 1);
    int total_seats = 0;
    for (int i = 1; i <= n; ++i) {
        cin >> parties[i].seats;
        parties[i].id = i;
        total_seats += parties[i].seats;
    }

    sort(parties.begin() + 1, parties.end(), [](const Party& a, const Party& b) {
        return a.seats > b.seats;
    });

    int target = total_seats / 2;

    dp[0][0] = true;
    int max_sum = -1;
    int best_i = -1;
    int best_s = -1;

    for (int i = 1; i <= n; ++i) {
        int w = parties[i].seats;

        for (int s = 0; s <= target; ++s) {
            if (dp[i - 1][s]) {
                if (s + w > target) {
                    if (s + w > max_sum) {
                        max_sum = s + w;
                        best_i = i;
                        best_s = s;
                    }
                }
            }
        }

        for (int s = 0; s <= target; ++s) {
            dp[i][s] = dp[i - 1][s];
            if (s >= w && dp[i - 1][s - w]) {
                dp[i][s] = true;
            }
        }
    }

    vector<int> result;
    if (max_sum != -1) {
        result.push_back(parties[best_i].id);
        int curr_s = best_s;

        for (int i = best_i - 1; i >= 1; --i) {
            if (curr_s == 0) break;
            if (curr_s >= parties[i].seats && dp[i - 1][curr_s - parties[i].seats]) {
                result.push_back(parties[i].id);
                curr_s -= parties[i].seats;
            }
        }

        sort(result.begin(), result.end());
    }

    cout << result.size() << "\n";
    for (size_t i = 0; i < result.size(); ++i) {
        cout << result[i] << (i + 1 == result.size() ? "" : " ");
    }
    cout << "\n";

    return 0;
}