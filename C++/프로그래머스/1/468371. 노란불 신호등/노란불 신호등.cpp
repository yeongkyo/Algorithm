#include <vector>

using namespace std;

long long gcd(long long a, long long b) {
    while (b != 0) {
        long long temp = a % b;
        a = b;
        b = temp;
    }
    return a;
}

long long lcm(long long a, long long b) {
    if (a == 0 || b == 0) return 0;
    return (a / gcd(a, b)) * b;
}

int solution(vector<vector<int>> signals) {
    long long L = 1;
    for (const auto& sig : signals) {
        long long T = sig[0] + sig[1] + sig[2];
        L = lcm(L, T);
    }

    for (long long t = 1; t <= L; ++t) {
        bool all_yellow = true;

        for (const auto& sig : signals) {
            int G = sig[0];
            int Y = sig[1];
            int R = sig[2];
            int T = G + Y + R;

            int rem = (t - 1) % T + 1;

            if (rem <= G || rem > G + Y) {
                all_yellow = false;
                break;
            }
        }

        if (all_yellow) {
            return t;
        }
    }

    return -1;
}