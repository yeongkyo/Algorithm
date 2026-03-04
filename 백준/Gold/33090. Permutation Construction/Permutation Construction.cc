#include <iostream>
#include <vector>
#include <stack>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int N;
    if (!(cin >> N)) return 0;

    vector<int> A(N + 1);
    vector<vector<int>> adj(N + 2); 

    for (int i = 1; i <= N; ++i) {
        cin >> A[i];
        if (A[i] == -1) {
            adj[N + 1].push_back(i);
        } else {
            adj[A[i]].push_back(i);
        }
    }

    vector<int> P(N + 1, 0);
    int current_val = N;
    
    stack<int> dfs_st;
    dfs_st.push(N + 1);

    while (!dfs_st.empty()) {
        int u = dfs_st.top();
        dfs_st.pop();

        // 루트 노드가 아닐 경우 순열 값 할당
        if (u != N + 1) {
            P[u] = current_val--;
        }

        // 작은 인덱스의 자식부터 먼저 방문해야 하므로, 스택에는 역순으로 삽입
        for (int i = (int)adj[u].size() - 1; i >= 0; --i) {
            dfs_st.push(adj[u][i]);
        }
    }

    vector<int> NGE(N + 1, -1);
    stack<int> mono_st;
    
    for (int i = 1; i <= N; ++i) {
        while (!mono_st.empty() && P[mono_st.top()] < P[i]) {
            NGE[mono_st.top()] = i;
            mono_st.pop();
        }
        mono_st.push(i);
    }

    bool possible = true;
    for (int i = 1; i <= N; ++i) {
        if (NGE[i] != A[i]) {
            possible = false;
            break;
        }
    }

    if (possible) {
        for (int i = 1; i <= N; ++i) {
            cout << P[i] << (i == N ? "" : " ");
        }
        cout << "\n";
    } else {
        cout << "-1\n";
    }

    return 0;
}