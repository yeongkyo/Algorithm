#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

struct Block {
    int digit;
    long long count;
};

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    if (!(cin >> n)) return 0;

    vector<long long> p(n);
    long long c_min = 2e18;
    int d_min = -1;

    for (int i = 0; i < n; i++) {
        cin >> p[i];
    }

    for (int i = n - 1; i >= 0; i--) {
        if (p[i] < c_min) {
            c_min = p[i];
            d_min = i;
        }
    }

    long long c_min_nz = 2e18;
    int d_min_nz = -1;
    for (int i = n - 1; i >= 1; i--) {
        if (p[i] < c_min_nz) {
            c_min_nz = p[i];
            d_min_nz = i;
        }
    }

    long long m;
    cin >> m;

    if (d_min_nz == -1 || m < c_min_nz) {
        if (m >= p[0]) {
            cout << "1\n0\n0\n";
        } else {
            cout << "0\n";
        }
        return 0;
    }

    long long L = 1 + (m - c_min_nz) / c_min;
    long long rem_L = L;
    vector<Block> blocks;

    while (rem_L > 0) {
        for (int d = n - 1; d >= 0; d--) {
            if (rem_L == L && d == 0) continue; 

            if (m >= p[d]) {
                long long rem_m = m - p[d];
                if (rem_m >= (rem_L - 1) * c_min) {
                    long long cnt;
                    if (p[d] == c_min) {
                        cnt = rem_L;
                    } else {
                        long long surplus = m - rem_L * c_min;
                        cnt = surplus / (p[d] - c_min);
                        if (cnt > rem_L) cnt = rem_L;
                    }

                    blocks.push_back({d, cnt});
                    m -= cnt * p[d];
                    rem_L -= cnt;
                    break;
                }
            }
        }
    }

    cout << L << "\n";

    if (L <= 50) {
        string s = "";
        for (auto& b : blocks) {
            s.append(b.count, '0' + b.digit);
        }
        cout << s << "\n" << s << "\n";
    } else {
        string first_50 = "";
        long long needed = 50;
        for (auto& b : blocks) {
            long long take = min(needed, b.count);
            first_50.append(take, '0' + b.digit);
            needed -= take;
            if (needed == 0) break;
        }

        string last_50 = "";
        needed = 50;
        for (int i = (int)blocks.size() - 1; i >= 0; i--) {
            long long take = min(needed, blocks[i].count);
            string temp(take, '0' + blocks[i].digit);
            last_50 = temp + last_50;
            needed -= take;
            if (needed == 0) break;
        }

        cout << first_50 << "\n" << last_50 << "\n";
    }

    return 0;
}