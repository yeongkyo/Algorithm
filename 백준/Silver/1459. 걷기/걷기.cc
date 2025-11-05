#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    
    long long X, Y, W, S;
    if (!(cin >> X >> Y >> W >> S)) return 0;
    
    long long a = min(X, Y);
    long long b = max(X, Y);
    long long ans = 0;
    
    if (S >= 2 * W) {
        ans = (X + Y) * W;
    } else if (S <= W) {
        if ((X + Y) % 2 == 0) ans = b * S;
        else ans = (b - 1) * S + W;
    } else {
        ans = a * S + (b - a) * W;
    }
    
    cout << ans << '\n';
    return 0;
}