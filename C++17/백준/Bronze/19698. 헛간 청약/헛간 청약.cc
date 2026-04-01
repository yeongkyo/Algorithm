#include <iostream>
#include <algorithm>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    long long N, W, H, L;
    cin >> N >> W >> H >> L;

    long long max_cows = (W / L) * (H / L);
    
    cout << min(N, max_cows) << "\n";

    return 0;
}