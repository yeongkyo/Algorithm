#include <iostream>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    cin >> n;

    int total_sum = 0;
    for (int i = 0; i < n; i++) {
        int move;
        cin >> move;
        total_sum += move;
    }

    if (total_sum > 0) {
        cout << "Right" << "\n";
    } else if (total_sum < 0) {
        cout << "Left" << "\n";
    } else {
        cout << "Stay" << "\n";
    }

    return 0;
}