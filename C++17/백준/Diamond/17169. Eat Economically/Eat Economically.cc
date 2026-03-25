#include <iostream>
#include <vector>
#include <algorithm>
#include <set>

using namespace std;

struct Item {
    long long l, d;
    int id;
    bool operator<(const Item& other) const {
        long long diff1 = l - d;
        long long diff2 = other.l - other.d;
        if (diff1 != diff2) return diff1 < diff2;
        return id < other.id;
    }
};

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    if (!(cin >> n)) return 0;

    vector<Item> items(2 * n);
    for (int i = 0; i < 2 * n; ++i) {
        cin >> items[i].l >> items[i].d;
        items[i].id = i;
    }

    sort(items.begin(), items.end());

    for (int i = 0; i < 2 * n; ++i) {
        items[i].id = i;
    }

    set<pair<long long, int>> L_vals; // {비용, 인덱스}
    set<int> left_active;             // 인덱스 모음
    set<pair<long long, int>> D_vals;
    set<int> right_active;

    long long current_cost = 0;
    
    for (int i = 0; i < n; ++i) {
        current_cost += items[i].l;
        left_active.insert(i);
        L_vals.insert({items[i].l, i});
    }
    for (int i = n; i < 2 * n; ++i) {
        current_cost += items[i].d;
        right_active.insert(i);
        D_vals.insert({items[i].d, i});
    }

    vector<long long> ans(n + 1);
    ans[n] = current_cost;

    for (int i = n; i >= 2; --i) {
        long long cost1 = -2e18, cost2 = -2e18, cost3 = -2e18;

        auto itL = L_vals.rbegin();
        pair<long long, int> pL1 = *itL; ++itL;
        pair<long long, int> pL2 = *itL;

        auto itD = D_vals.rbegin();
        pair<long long, int> pD1 = *itD; ++itD;
        pair<long long, int> pD2 = *itD;

        int first_R = *right_active.begin();
        int last_L = *left_active.rbegin();

        cost1 = pL1.first + pD1.first;
        cost2 = pL1.first + pL2.first + items[first_R].d - items[first_R].l;
        cost3 = pD1.first + pD2.first + items[last_L].l - items[last_L].d;

        long long best = max({cost1, cost2, cost3});
        current_cost -= best;
        ans[i - 1] = current_cost;

        if (best == cost1) {
            L_vals.erase(pL1);
            left_active.erase(pL1.second);
            D_vals.erase(pD1);
            right_active.erase(pD1.second);
        } else if (best == cost2) {
            L_vals.erase(pL1);
            left_active.erase(pL1.second);
            L_vals.erase(pL2);
            left_active.erase(pL2.second);

            right_active.erase(first_R);
            D_vals.erase({items[first_R].d, first_R});

            left_active.insert(first_R);
            L_vals.insert({items[first_R].l, first_R});
        } else {
            D_vals.erase(pD1);
            right_active.erase(pD1.second);
            D_vals.erase(pD2);
            right_active.erase(pD2.second);

            left_active.erase(last_L);
            L_vals.erase({items[last_L].l, last_L});

            right_active.insert(last_L);
            D_vals.insert({items[last_L].d, last_L});
        }
    }

    for (int i = 1; i <= n; ++i) {
        cout << ans[i] << "\n";
    }

    return 0;
}