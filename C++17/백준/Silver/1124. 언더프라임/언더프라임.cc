#include <iostream>
#include <vector>

using namespace std;

const int MAX = 100000;
int spf[MAX + 1];
int factors_count[MAX + 1];

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int A, B;
    if (!(cin >> A >> B)) return 0;

    for (int i = 2; i <= MAX; i++) {
        spf[i] = i;
    }

    for (int i = 2; i * i <= MAX; i++) {
        if (spf[i] == i) {
            for (int j = i * i; j <= MAX; j += i) {
                if (spf[j] == j) {
                    spf[j] = i;
                }
            }
        }
    }

    for (int i = 2; i <= MAX; i++) {
        factors_count[i] = factors_count[i / spf[i]] + 1;
    }

    int ans = 0;
    for (int i = A; i <= B; i++) {
        int c = factors_count[i];
        if (c >= 2 && spf[c] == c) {
            ans++;
        }
    }

    cout << ans << "\n";

    return 0;
}