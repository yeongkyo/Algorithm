#include <iostream>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    cin >> n;

    int even_count = 0;
    int odd_count = 0;

    for (int i = 0; i < n; i++) {
        int a;
        cin >> a;
        if (a % 2 == 0) {
            even_count++;
        } else {
            odd_count++;
        }
    }

    if (even_count > odd_count) {
        cout << "Happy" << "\n";
    } else {
        cout << "Sad" << "\n";
    }

    return 0;
}