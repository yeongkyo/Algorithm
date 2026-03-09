#include <iostream>
#include <string>
#include <cmath>

using namespace std;

long long gcd(long long a, long long b) {
    return b == 0 ? a : gcd(b, a % b);
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    long long a, b, c, d;
    char slash;

    if (!(cin >> a >> slash >> b)) return 0;
    if (!(cin >> c >> slash >> d)) return 0;

    long long num = a * d - c * b;
    long long den = b * d;

    long long common = gcd(abs(num), den);
    num /= common;
    den /= common;

    long long int_part = num / den;
    long long frac_num = abs(num % den);
    long long frac_den = den;

    if (num < 0 && int_part == 0) {
        cout << "-0\n";
    } else {
        cout << int_part << "\n";
    }

    if (frac_num != 0) {
        cout << frac_num << "/" << frac_den << "\n";

        string int_str;
        if (num < 0 && int_part == 0) {
            int_str = "-";
        } else if (int_part == 0) {
            int_str = "";
        } else {
            int_str = to_string(int_part);
        }

        string s_num = to_string(frac_num);
        string s_den = to_string(frac_den);
        int dash_len = (int)s_den.length();
        int int_len = (int)int_str.length();

        for (int i = 0; i < int_len + dash_len - (int)s_num.length(); ++i) {
            cout << " ";
        }
        cout << s_num << "\n";

        cout << int_str;
        for (int i = 0; i < dash_len; ++i) {
            cout << "-";
        }
        cout << "\n";

        for (int i = 0; i < int_len; ++i) {
            cout << " ";
        }
        cout << s_den << "\n";
    } else {
        cout << "\n";
    }

    return 0;
}