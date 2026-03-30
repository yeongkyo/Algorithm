#include <iostream>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    long long N, M;
    if (!(cin >> N >> M)) return 0;

    for (int i = 0; i < N; ++i) {
        long long a;
        cin >> a;
    }

    long long phi = M;
    long long temp = M;
    
    for (long long i = 2; i * i <= temp; ++i) {
        if (temp % i == 0) {
            while (temp % i == 0) {
                temp /= i;
            }
            phi -= phi / i;
        }
    }
    
    if (temp > 1) {
        phi -= phi / temp;
    }

    cout << phi << "\n";
    return 0;
}