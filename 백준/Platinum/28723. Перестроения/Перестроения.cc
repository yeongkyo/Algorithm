#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;


int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    if (!(cin >> n)) return 0;

    vector<int> a(n), b(n);
    for (int i = 0; i < n; ++i) cin >> a[i];
    for (int i = 0; i < n; ++i) cin >> b[i];

    // 현재 상태와 목표 상태가 같으면 연산 필요 없음
    if (a == b) {
        cout << 0 << endl;
        return 0;
    }

    vector<int> target_pos(n + 1);
    for (int i = 0; i < n; ++i) {
        target_pos[b[i]] = i;
    }

    vector<vector<int>> results;
    int max_bit = 0;
    while ((1 << max_bit) < n) max_bit++;

    for (int k = 0; k < max_bit; ++k) {
        vector<int> S;
        vector<bool> in_S(n + 1, false);
        
        for (int i = 0; i < n; ++i) {
            int v = a[i];
            int pos = target_pos[v];
            if (((pos >> k) & 1) != ((pos >> (k + 1)) & 1)) {
                S.push_back(v);
                in_S[v] = true;
            }
        }

        if (!S.empty()) {
            results.push_back(S);
            
            vector<int> subset_order;
            for (int i = 0; i < n; ++i) {
                if (in_S[a[i]]) subset_order.push_back(a[i]);
            }
            reverse(subset_order.begin(), subset_order.end());

            vector<int> next_a;
            next_a.reserve(n);
            for (int x : subset_order) next_a.push_back(x);
            for (int i = 0; i < n; ++i) {
                if (!in_S[a[i]]) next_a.push_back(a[i]);
            }
            a = next_a;
        }
    }

    bool already_done = true;
    for (int i = 0; i < n; ++i) if (a[i] != b[i]) already_done = false;
    
    if (!already_done) {
        results.push_back(a);
    }

    cout << results.size() << "\n";
    for (const auto& s_vec : results) {
        cout << s_vec.size();
        for (int friend_id : s_vec) {
            cout << " " << friend_id;
        }
        cout << "\n";
    }

    return 0;
}