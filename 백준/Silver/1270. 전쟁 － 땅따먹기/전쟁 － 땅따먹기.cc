#include <iostream>
#include <map>

using namespace std;

void solve() {
    long long t;
    if (!(cin >> t)) return;

    map<long long, long long> counts;
    long long winner = -1;
    bool found = false;

    for (int i = 0; i < t; i++) {
        long long army_id;
        cin >> army_id;
        
        counts[army_id]++;
        
        if (counts[army_id] > t / 2) {
            if (!found) {
                winner = army_id;
                found = true;
            }
        }
    }

    if (found) {
        cout << winner << "\n";
    } else {
        cout << "SYJKGW\n";
    }
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    if (!(cin >> n)) return 0;

    while (n--) {
        solve();
    }

    return 0;
}