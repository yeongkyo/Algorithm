#include <bits/stdc++.h>
using namespace std;

int W, H;

int convert(int dir, int dist) {
    if (dir == 1) return dist;
    if (dir == 4) return W + dist;
    if (dir == 2) return W + H + (W - dist);
    return W + H + W + (H - dist);
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    cin >> W >> H;

    int N;
    cin >> N;

    vector<int> pos(N);
    for (int i = 0; i < N; i++) {
        int d, x;
        cin >> d >> x;
        pos[i] = convert(d, x);
    }

    int d, x;
    cin >> d >> x;
    int me = convert(d, x);

    int perimeter = 2 * (W + H);
    int ans = 0;

    for (int p : pos) {
        int diff = abs(me - p);
        ans += min(diff, perimeter - diff);
    }

    cout << ans << '\n';
    return 0;
}