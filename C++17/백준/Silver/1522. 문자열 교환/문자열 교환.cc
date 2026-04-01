#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    string s;
    cin >> s;

    int n = s.length();
    int countA = 0;
    for (char c : s) {
        if (c == 'a') countA++;
    }

    if (countA == 0 || countA == n) {
        cout << 0 << "\n";
        return 0;
    }

    string extended = s + s;
    int currentB = 0;
    for (int i = 0; i < countA; i++) {
        if (extended[i] == 'b') currentB++;
    }

    int minSwaps = currentB;

    for (int i = 1; i < n; i++) {
        if (extended[i - 1] == 'b') currentB--;
        if (extended[i + countA - 1] == 'b') currentB++;
        minSwaps = min(minSwaps, currentB);
    }

    cout << minSwaps << "\n";

    return 0;
}