#include <iostream>
#include <vector>

using namespace std;

const long long INF = 2e18;

struct Line {
    long long a, b;
    long long get(long long x) {
        return a * x + b;
    }
};

struct Node {
    int left, right;
    Line line;
};

vector<Node> tree;

void init() {
    tree.push_back({-1, -1, {0, INF}});
}

void insert(int node, long long start, long long end, Line newLine) {
    long long mid = (start + end) / 2;
    Line low = tree[node].line;
    Line high = newLine;
    
    if (low.get(start) > high.get(start)) swap(low, high);
    
    if (low.get(end) <= high.get(end)) {
        tree[node].line = low;
        return;
    }
    
    if (low.get(mid) < high.get(mid)) {
        tree[node].line = low;
        if (tree[node].right == -1) {
            tree[node].right = tree.size();
            tree.push_back({-1, -1, {0, INF}});
        }
        insert(tree[node].right, mid + 1, end, high);
    } else {
        tree[node].line = high;
        if (tree[node].left == -1) {
            tree[node].left = tree.size();
            tree.push_back({-1, -1, {0, INF}});
        }
        insert(tree[node].left, start, mid, low);
    }
}

long long query(int node, long long start, long long end, long long x) {
    if (node == -1) return INF;
    long long mid = (start + end) / 2;
    long long res = tree[node].line.get(x);
    if (x <= mid) {
        res = min(res, query(tree[node].left, start, mid, x));
    } else {
        res = min(res, query(tree[node].right, mid + 1, end, x));
    }
    return res;
}

int n;
long long dist[100005];
long long p[100005], s[100005];
long long dp[100005];

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    if (!(cin >> n)) return 0;

    for (int i = 2; i <= n; i++) {
        long long d;
        cin >> d;
        dist[i] = dist[i - 1] + d;
    }

    for (int i = 1; i < n; i++) {
        cin >> p[i] >> s[i];
    }

    init();
    dp[1] = 0;
    insert(0, 0, dist[n], {s[1], dp[1] + p[1] - s[1] * dist[1]});

    for (int i = 2; i <= n; i++) {
        dp[i] = query(0, 0, dist[n], dist[i]);
        if (i < n) {
            insert(0, 0, dist[n], {s[i], dp[i] + p[i] - s[i] * dist[i]});
        }
    }

    cout << dp[n] << "\n";

    return 0;
}