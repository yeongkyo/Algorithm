#include <bits/stdc++.h>
using namespace std;

struct Compare {
    bool operator()(int a, int b) const {
        if (abs(a) == abs(b)) return a > b;
        return abs(a) > abs(b);
    }
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int N;
    cin >> N;

    priority_queue<int, vector<int>, Compare> pq;

    while (N--) {
        int x;
        cin >> x;

        if (x != 0) {
            pq.push(x);
        } else {
            if (pq.empty()) {
                cout << 0 << '\n';
            } else {
                cout << pq.top() << '\n';
                pq.pop();
            }
        }
    }

    return 0;
}