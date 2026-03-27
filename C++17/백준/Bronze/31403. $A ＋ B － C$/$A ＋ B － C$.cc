#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int A, B, C;
    cin >> A >> B >> C;

    cout << A + B - C << '\n';

    string sA = to_string(A);
    string sB = to_string(B);
    string sC = to_string(C);

    string concat = sA + sB;
    cout << stoi(concat) - C << '\n';

    return 0;
}