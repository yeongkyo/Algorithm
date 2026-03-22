#include <iostream>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    cin >> n;

    int total_months = 8 + (n - 1) * 7;
    int year = 2024 + (total_months - 1) / 12;
    int month = (total_months - 1) % 12 + 1;

    cout << year << " " << month << "\n";

    return 0;
}