#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

struct Item {
    long long w, v;
};

int N;
long long C;
vector<Item> items;
vector<long long> rem_w, rem_v;
long long ans = 0;

void dfs(int idx, long long cw, long long cv, bool prev_skipped) {
    if (cv > ans) ans = cv;
    if (idx == N) return;
    
    if (cv + rem_v[idx] <= ans) return;
    if (cw + rem_w[idx] <= C) {
        if (cv + rem_v[idx] > ans) {
            ans = cv + rem_v[idx];
        }
        return;
    }

    if (cw + items[idx].w <= C) {
        bool dominated = false;
        if (idx > 0 && items[idx].w == items[idx - 1].w && prev_skipped) {
            dominated = true;
        }
        if (!dominated) {
            dfs(idx + 1, cw + items[idx].w, cv + items[idx].v, false);
        }
    }

    dfs(idx + 1, cw, cv, true);
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    if (!(cin >> N)) return 0;
    
    items.resize(N);
    for (int i = 0; i < N; ++i) {
        cin >> items[i].w >> items[i].v;
    }
    cin >> C;
    
    sort(items.begin(), items.end(), [](const Item& a, const Item& b) {
        if (a.w != b.w) return a.w > b.w;
        return a.v > b.v;
    });
    
    rem_w.assign(N + 1, 0);
    rem_v.assign(N + 1, 0);
    for (int i = N - 1; i >= 0; --i) {
        rem_w[i] = rem_w[i + 1] + items[i].w;
        rem_v[i] = rem_v[i + 1] + items[i].v;
    }
    
    dfs(0, 0, 0, false);
    
    cout << ans << "\n";
    return 0;
}