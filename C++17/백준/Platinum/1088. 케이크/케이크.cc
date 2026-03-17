#include <iostream>
#include <vector>
#include <queue>
#include <iomanip>

using namespace std;

struct Piece {
    long long w;
    long long c;
    int id;
    
    bool operator<(const Piece& other) const {
        return w * other.c < other.w * c;
    }
};

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    if (!(cin >> n)) return 0;

    vector<long long> w(n);
    long long w_min = -1;
    long long c_min = 1;
    
    priority_queue<Piece> pq;

    for (int i = 0; i < n; ++i) {
        cin >> w[i];
        pq.push({w[i], 1, i});
        if (w_min == -1 || w[i] < w_min) {
            w_min = w[i];
            c_min = 1;
        }
    }

    int m;
    cin >> m;

    double ans = (double)pq.top().w / pq.top().c - (double)w_min / c_min;

    for (int i = 0; i < m; ++i) {
        Piece p = pq.top();
        pq.pop();

        p.c++;
        pq.push(p);

        if (p.w * c_min < w_min * p.c) {
            w_min = p.w;
            c_min = p.c;
        }

        double cur_max = (double)pq.top().w / pq.top().c;
        double cur_min = (double)w_min / c_min;
        
        if (cur_max - cur_min < ans) {
            ans = cur_max - cur_min;
        }
    }

    cout << fixed << setprecision(10) << ans << "\n";
    return 0;
}