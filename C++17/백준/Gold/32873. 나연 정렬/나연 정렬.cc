#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int N;
    if (!(cin >> N)) return 0;

    vector<int> tails;
    
    for (int i = 0; i < N; i++) {
        int x;
        cin >> x;
        
        if (tails.empty() || x > tails.back()) {
            tails.push_back(x);
        } else {
            auto it = lower_bound(tails.begin(), tails.end(), x);
            *it = x;
        }
    }

    cout << tails.size() << "\n";

    return 0;
}