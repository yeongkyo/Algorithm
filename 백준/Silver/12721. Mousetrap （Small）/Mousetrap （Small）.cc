#include <iostream>
#include <vector>

using namespace std;

void solve(int tc) {
    int K;
    cin >> K;
    
    vector<int> empty_spots(K);
    for (int i = 0; i < K; ++i) {
        empty_spots[i] = i + 1;
    }
    
    vector<int> deck(K + 1, 0);
    int curr = 0; 
    
    for (int i = 1; i <= K; ++i) {
        curr = (curr + i - 1) % empty_spots.size();
        
        deck[empty_spots[curr]] = i;
        
        empty_spots.erase(empty_spots.begin() + curr);
    }
    
    int n;
    cin >> n;
    cout << "Case #" << tc << ":";
    for (int i = 0; i < n; ++i) {
        int d;
        cin >> d;
        cout << " " << deck[d];
    }
    cout << "\n";
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    int T;
    if (cin >> T) {
        for (int i = 1; i <= T; ++i) {
            solve(i);
        }
    }
    return 0;
}