#include <iostream>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    long long n, k;
    if (!(cin >> n >> k)) return 0;

    if (n % 2 != 0 || k > n || k == 1) {
        cout << -1 << "\n";
        return 0;
    }

    if (k % 2 == 0) {
        for (int i = 0; i < k - 2; ++i) {
            cout << 1 << " ";
        }
        long long x = (n - k + 2) / 2;
        cout << x << " " << x << "\n";
    } else {
        if (n == k + 1) {
            cout << -1 << "\n";
            return 0;
        }
        
        long long S = n - k + 3;
        bool is_power_of_2 = (S & (S - 1)) == 0;

        if (is_power_of_2) {
            if (k == 3) {
                cout << -1 << "\n";
                return 0;
            }
            
            for (int i = 0; i < k - 5; ++i) {
                cout << 1 << " ";
            }
            cout << 2 << " " << 2 << " ";
            
            long long X = S - 2;
            long long half = X / 2;
            long long A = half & -half;
            long long B = half - A;
            long long C = half;
            
            cout << A << " " << B << " " << C << "\n";
            
        } else {
            for (int i = 0; i < k - 3; ++i) {
                cout << 1 << " ";
            }
            
            long long half = S / 2;
            long long A = half & -half;
            long long B = half - A;
            long long C = half;
            
            cout << A << " " << B << " " << C << "\n";
        }
    }

    return 0;
}