#include <iostream>
#include <string>
#include <algorithm>

using namespace std;

// 숫자를 뒤집어서 반환하는 함수
long long get_reversed(long long n) {
    long long reversed_n = 0;
    while (n > 0) {
        reversed_n = reversed_n * 10 + (n % 10);
        n /= 10;
    }
    return reversed_n;
}

bool is_palindrome(long long n) {
    string s = to_string(n);
    string original = s;
    reverse(s.begin(), s.end()); 
    return original == s;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    int T;
    cin >> T;

    while (T--) {
        long long N;
        cin >> N;

        long long sum_val = N + get_reversed(N);

        if (is_palindrome(sum_val)) {
            cout << "YES\n";
        } else {
            cout << "NO\n";
        }
    }

    return 0;
}