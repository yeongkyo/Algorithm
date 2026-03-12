#include <iostream>
#include <string>
#include <vector>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    string A, B;
    if (!(cin >> A >> B)) return 0;

    int n = A.length();
    int m = B.length();

    vector<int> V(m + 1);
    for (int k = 1; k <= m; k++) {
        V[k] = k;
    }

    long long total_sum = 0;

    for (int i = 0; i < n; i++) {
        int L = 0;
        long long current_S = 0; 
        
        for (int k = 1; k <= m; k++) {
            int temp = V[k];
            
            if (A[i] == B[k - 1]) {
                V[k] = L;
                L = temp;
            } 
            else if (V[k] < L) {
                V[k] = L;
                L = temp;
            }
            
            current_S += (k - V[k]);
            
            total_sum += current_S;
        }
    }

    cout << total_sum << "\n";

    return 0;
}