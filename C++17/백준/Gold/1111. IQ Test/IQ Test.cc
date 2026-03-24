#include <iostream>
#include <vector>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int N;
    if (!(cin >> N)) return 0;

    vector<int> arr(N);
    for (int i = 0; i < N; ++i) {
        cin >> arr[i];
    }

    if (N == 1) {
        cout << "A\n";
    } else if (N == 2) {
        if (arr[0] == arr[1]) {
            cout << arr[0] << "\n";
        } else {
            cout << "A\n";
        }
    } else {
        int a = 0, b = 0;
        if (arr[0] == arr[1]) {
            a = 1;
            b = 0;
        } else {
            if ((arr[2] - arr[1]) % (arr[1] - arr[0]) != 0) {
                cout << "B\n";
                return 0;
            }
            a = (arr[2] - arr[1]) / (arr[1] - arr[0]);
            b = arr[1] - arr[0] * a;
        }

        bool possible = true;
        for (int i = 1; i < N; ++i) {
            if (arr[i] != arr[i - 1] * a + b) {
                possible = false;
                break;
            }
        }

        if (possible) {
            cout << arr[N - 1] * a + b << "\n";
        } else {
            cout << "B\n";
        }
    }

    return 0;
}