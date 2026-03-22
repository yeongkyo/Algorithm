#include <iostream>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    long long t1, e1, f1;
    long long t2, e2, f2;

    cin >> t1 >> e1 >> f1;
    cin >> t2 >> e2 >> f2;

    long long max_time = t1 * 3 + e1 * 20 + f1 * 120;
    long long mel_time = t2 * 3 + e2 * 20 + f2 * 120;

    if (max_time > mel_time) {
        cout << "Max" << "\n";
    } else if (mel_time > max_time) {
        cout << "Mel" << "\n";
    } else {
        cout << "Draw" << "\n";
    }

    return 0;
}