#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

typedef __int128_t int128;

int n;
long long dp_cnt[20][120];
int128 dp_sum[20][120];
int128 p10[20];
long long p10_n = 1;

void precompute() {
    p10[0] = 1;
    for (int i = 1; i <= 18; i++) {
        p10[i] = p10[i - 1] * 10;
    }
    for (int i = 0; i < n; i++) {
        p10_n *= 10;
    }

    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= 108; j++) {
            dp_cnt[i][j] = 0;
            dp_sum[i][j] = 0;
        }
    }

    dp_cnt[n][54] = 1;

    for (int pos = n - 1; pos >= 0; pos--) {
        for (int diff = 0; diff <= 108; diff++) {
            for (int d = 0; d <= 9; d++) {
                int nxt_diff = diff + (pos < n / 2 ? d : -d);
                if (nxt_diff >= 0 && nxt_diff <= 108) {
                    long long cnt = dp_cnt[pos + 1][nxt_diff];
                    int128 sval = dp_sum[pos + 1][nxt_diff];
                    if (cnt > 0) {
                        dp_cnt[pos][diff] += cnt;
                        dp_sum[pos][diff] += sval + (int128)cnt * d * p10[n - 1 - pos];
                    }
                }
            }
        }
    }
}

pair<long long, int128> get_f_h(string R_str) {
    long long ans_cnt = 0;
    int128 ans_sum = 0;
    int curr_diff = 54;
    long long prefix_val = 0;

    for (int pos = 0; pos < n; pos++) {
        int limit = R_str[pos] - '0';
        for (int d = 0; d < limit; d++) {
            int nxt_diff = curr_diff + (pos < n / 2 ? d : -d);
            if (nxt_diff >= 0 && nxt_diff <= 108) {
                long long cnt = dp_cnt[pos + 1][nxt_diff];
                int128 sval = dp_sum[pos + 1][nxt_diff];
                if (cnt > 0) {
                    ans_cnt += cnt;
                    ans_sum += sval + (int128)cnt * (prefix_val * 10 + d) * p10[n - 1 - pos];
                }
            }
        }
        curr_diff += (pos < n / 2 ? limit : -limit);
        prefix_val = prefix_val * 10 + limit;
    }

    if (curr_diff == 54) {
        ans_cnt += 1;
        ans_sum += prefix_val;
    }

    return {ans_cnt, ans_sum};
}

string to_ndigit_string(long long val) {
    string s = to_string(val);
    while (s.length() < n) {
        s = "0" + s;
    }
    return s;
}

long long get_W(long long R) {
    if (R < 0) return 0;
    return get_f_h(to_ndigit_string(R)).first;
}

int128 get_g(long long R) {
    if (R < 0) return 0;
    auto result = get_f_h(to_ndigit_string(R));
    long long f = result.first;
    int128 h = result.second;
    return (int128)(R + 1) * f - h;
}

int128 get_G(long long X) {
    if (X < 0) return 0;
    long long Q = X / p10_n;
    long long R = X % p10_n;

    int128 Q_128 = Q;
    int128 C = get_W(p10_n - 1);
    int128 S_W = get_g(p10_n - 1);

    int128 ans = p10[n] * C * Q_128 * (Q_128 - 1) / 2;
    ans += Q_128 * S_W;
    ans += (R + 1) * Q_128 * C;
    ans += get_g(R);
    return ans;
}

long long gcd(long long a, long long b) {
    while (b) {
        a %= b;
        swap(a, b);
    }
    return a;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    string a_str, b_str;
    long long k;

    if (cin >> a_str >> b_str >> k) {
        n = a_str.length(); 
        precompute();

        long long a = stoll(a_str);
        long long b = stoll(b_str);

        int128 total_lucky_sum = get_G(b + k - 1) - get_G(b - 1) - get_G(a + k - 2) + get_G(a - 2);
        
        long long numerator = (long long)total_lucky_sum;
        long long denominator = b - a + 1;

        long long g = gcd(numerator, denominator);
        numerator /= g;
        denominator /= g;

        if (denominator == 1) {
            cout << numerator << "\n";
        } else {
            cout << numerator << "/" << denominator << "\n";
        }
    }

    return 0;
}