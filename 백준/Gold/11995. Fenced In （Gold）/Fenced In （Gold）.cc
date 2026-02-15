#include <bits/stdc++.h>
using namespace std;

using ll = long long;

struct DSU {
    int n;
    vector<int> p, r;
    DSU(int n=0): n(n), p(n), r(n,0) {
        iota(p.begin(), p.end(), 0);
    }
    int find(int x){
        while(p[x]!=x){
            p[x]=p[p[x]];
            x=p[x];
        }
        return x;
    }
    bool unite(int a,int b){
        a=find(a); b=find(b);
        if(a==b) return false;
        if(r[a]<r[b]) swap(a,b);
        p[b]=a;
        if(r[a]==r[b]) r[a]++;
        return true;
    }
};

struct Item {
    int type; // 0: column(dx), 1: row(dy)
    int idx;
    int w;
    bool operator<(Item const& other) const {
        return w < other.w;
    }
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int A,B,n,m;
    cin >> A >> B >> n >> m;

    vector<int> xs(n+2), ys(m+2);
    xs[0]=0; xs[n+1]=A;
    ys[0]=0; ys[m+1]=B;
    for(int i=1;i<=n;i++) cin >> xs[i];
    for(int i=1;i<=m;i++) cin >> ys[i];
    sort(xs.begin(), xs.end());
    sort(ys.begin(), ys.end());

    // dx: column widths (n+1), dy: row heights (m+1)
    vector<int> dx(n+1), dy(m+1);
    for(int i=0;i<n+1;i++) dx[i]=xs[i+1]-xs[i];
    for(int i=0;i<m+1;i++) dy[i]=ys[i+1]-ys[i];

    // items: each column and each row
    vector<Item> items;
    items.reserve((n+1)+(m+1));
    for(int i=0;i<n+1;i++) items.push_back({0,i,dx[i]});
    for(int j=0;j<m+1;j++) items.push_back({1,j,dy[j]});
    sort(items.begin(), items.end());

    int cols = n+1;
    int rows = m+1;
    int V = cols * rows;
    auto id = [&](int c,int r){ return c*rows + r; };

    DSU dsu(V);
    ll ans = 0;
    int added = 0;
    int need = V - 1;

    for (auto &it : items) {
        if (added == need) break;

        if (it.type == 0) {
            // column it.idx: connect (r) <-> (r+1) for r=0..m-1, cost = dx[col]
            int c = it.idx;
            for (int r = 0; r < m; r++) {
                if (dsu.unite(id(c,r), id(c,r+1))) {
                    ans += (ll)it.w;
                    if (++added == need) break;
                }
            }
        } else {
            // row it.idx: connect (c) <-> (c+1) for c=0..n-1, cost = dy[row]
            int r = it.idx;
            for (int c = 0; c < n; c++) {
                if (dsu.unite(id(c,r), id(c+1,r))) {
                    ans += (ll)it.w;
                    if (++added == need) break;
                }
            }
        }
    }

    cout << ans << "\n";
    return 0;
}
