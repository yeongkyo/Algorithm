#include <iostream>
#include <vector>
#include <algorithm>
#include <set>

using namespace std;

const int MAXN = 200005;
const int MAXC = 400010;

int S[MAXN], E[MAXN];
int jump[19][MAXC];
vector<int> coords;

int get_MIS(int L, int R) {
    if (L > R) return 0;
    int curr = L;
    int count = 0;
    for (int k = 18; k >= 0; --k) {
        if (curr <= R && jump[k][curr] <= R) {
            curr = jump[k][curr] + 1;
            count += (1 << k);
        }
    }
    return count;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int N;
    if (!(cin >> N)) return 0;

    coords.reserve(N * 2);
    for (int i = 1; i <= N; ++i) {
        cin >> S[i] >> E[i];
        coords.push_back(S[i]);
        coords.push_back(E[i]);
    }

    sort(coords.begin(), coords.end());
    coords.erase(unique(coords.begin(), coords.end()), coords.end());

    int C = coords.size();
    for (int i = 1; i <= N; ++i) {
        S[i] = lower_bound(coords.begin(), coords.end(), S[i]) - coords.begin() + 1;
        E[i] = lower_bound(coords.begin(), coords.end(), E[i]) - coords.begin() + 1;
    }

    for (int i = 1; i <= C + 2; ++i) jump[0][i] = C + 1;
    
    for (int i = 1; i <= N; ++i) {
        if (E[i] < jump[0][S[i]]) {
            jump[0][S[i]] = E[i];
        }
    }

    for (int i = C; i >= 1; --i) {
        if (jump[0][i + 1] < jump[0][i]) {
            jump[0][i] = jump[0][i + 1];
        }
    }

    for (int k = 1; k <= 18; ++k) {
        for (int i = 1; i <= C + 2; ++i) {
            int nxt = jump[k - 1][i] + 1;
            if (nxt > C + 1) nxt = C + 1;
            jump[k][i] = jump[k - 1][nxt];
        }
    }

    set<pair<int, int>> gaps;
    gaps.insert({1, C});

    vector<int> ans;
    ans.reserve(N);

    for (int i = 1; i <= N; ++i) {
        auto it = gaps.upper_bound({S[i], 1000000000});
        if (it == gaps.begin()) continue;
        --it;
        int L = it->first;
        int R = it->second;

        if (L <= S[i] && E[i] <= R) {
            int mis_all = get_MIS(L, R);
            int mis_left = get_MIS(L, S[i] - 1);
            int mis_right = get_MIS(E[i] + 1, R);

            if (mis_all == mis_left + 1 + mis_right) {
                ans.push_back(i);
                gaps.erase(it);
                if (L <= S[i] - 1) gaps.insert({L, S[i] - 1});
                if (E[i] + 1 <= R) gaps.insert({E[i] + 1, R});
            }
        }
    }

    cout << ans.size() << "\n";
    for (int i = 0; i < ans.size(); ++i) {
        cout << ans[i] << (i + 1 == ans.size() ? "" : " ");
    }
    cout << "\n";

    return 0;
}