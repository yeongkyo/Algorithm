#include <iostream>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int r, g, b;
    cin >> r >> g >> b;

    int total_cost = (r * 3) + (g * 4) + (b * 5);

    cout << total_cost << "\n";

    return 0;
}