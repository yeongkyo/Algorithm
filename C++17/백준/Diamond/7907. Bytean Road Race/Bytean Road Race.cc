#include <iostream>

using namespace std;

const int MAXN = 100005;
int X[MAXN], Y[MAXN];
int E_out[MAXN], S_out[MAXN];
int post1[MAXN], post2[MAXN];
bool visited1[MAXN], visited2[MAXN];
int timer1 = 0, timer2 = 0;

void dfs1(int u) {
    visited1[u] = true;
    if (S_out[u] && !visited1[S_out[u]]) dfs1(S_out[u]);
    if (E_out[u] && !visited1[E_out[u]]) dfs1(E_out[u]);
    post1[u] = ++timer1;
}

void dfs2(int u) {
    visited2[u] = true;
    if (E_out[u] && !visited2[E_out[u]]) dfs2(E_out[u]);
    if (S_out[u] && !visited2[S_out[u]]) dfs2(S_out[u]);
    post2[u] = ++timer2;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n, m, k;
    if (!(cin >> n >> m >> k)) return 0;

    for (int i = 1; i <= n; ++i) {
        cin >> X[i] >> Y[i];
    }

    for (int i = 0; i < m; ++i) {
        int u, v;
        cin >> u >> v;
        // 남쪽(South) 간선 찾기: X좌표가 같고 Y좌표가 변하는 경우
        if (X[u] == X[v]) {
            if (Y[u] > Y[v]) S_out[u] = v;
            else S_out[v] = u;
        } 
        // 동쪽(East) 간선 찾기: Y좌표가 같고 X좌표가 변하는 경우
        else {
            if (X[u] < X[v]) E_out[u] = v;
            else E_out[v] = u;
        }
    }

    dfs1(1);
    dfs2(1);

    for (int i = 0; i < k; ++i) {
        int p, q;
        cin >> p >> q;
        
        bool p_reaches_q = (post1[p] > post1[q]) && (post2[p] > post2[q]);
        bool q_reaches_p = (post1[q] > post1[p]) && (post2[q] > post2[p]);

        if (p_reaches_q || q_reaches_p) {
            cout << "TAK\n";
        } else {
            cout << "NIE\n";
        }
    }

    return 0;
}