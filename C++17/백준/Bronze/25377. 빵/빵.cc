#include <iostream>
#include <algorithm>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int N;
    cin >> N;

    int min_time = 2001;
    bool possible = false;

    for (int i = 0; i < N; i++) {
        int A, B;
        cin >> A >> B;

        if (A <= B) {
            possible = true;
            min_time = min(min_time, B);
        }
    }

    if (possible) {
        cout << min_time << "\n";
    } else {
        cout << -1 << "\n";
    }

    return 0;
}