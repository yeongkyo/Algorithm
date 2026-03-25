#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

struct Node {
    int val, cnt;
};

int n, c, m;
int arr[300005];
vector<int> pos[10005];
Node tree[1200005];

Node merge(Node a, Node b) {
    if (a.cnt == 0) return b;
    if (b.cnt == 0) return a;
    if (a.val == b.val) return {a.val, a.cnt + b.cnt};
    if (a.cnt > b.cnt) return {a.val, a.cnt - b.cnt};
    return {b.val, b.cnt - a.cnt};
}

void init(int node, int start, int end) {
    if (start == end) {
        tree[node] = {arr[start], 1};
        return;
    }
    int mid = (start + end) / 2;
    init(node * 2, start, mid);
    init(node * 2 + 1, mid + 1, end);
    tree[node] = merge(tree[node * 2], tree[node * 2 + 1]);
}

Node query(int node, int start, int end, int left, int right) {
    if (left > end || right < start) return {0, 0};
    if (left <= start && end <= right) return tree[node];
    int mid = (start + end) / 2;
    return merge(query(node * 2, start, mid, left, right),
                 query(node * 2 + 1, mid + 1, end, left, right));
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    if (!(cin >> n >> c)) return 0;
    for (int i = 1; i <= n; i++) {
        cin >> arr[i];
        pos[arr[i]].push_back(i);
    }

    init(1, 1, n);

    cin >> m;
    while (m--) {
        int a, b;
        cin >> a >> b;
        Node res = query(1, 1, n, a, b);
        int cand = res.val;
        int count = 0;

        if (cand > 0) {
            count = upper_bound(pos[cand].begin(), pos[cand].end(), b) -
                    lower_bound(pos[cand].begin(), pos[cand].end(), a);
        }

        if (count > (b - a + 1) / 2) {
            cout << "yes " << cand << "\n";
        } else {
            cout << "no\n";
        }
    }

    return 0;
}