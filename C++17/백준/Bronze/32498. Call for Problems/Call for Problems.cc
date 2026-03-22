#include <iostream>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    if (!(cin >> n)) return 0;

    int excluded_count = 0;
    for (int i = 0; i < n; ++i) {
        int d;
        cin >> d;
        if (d % 2 != 0) {
            excluded_count++;
        }
    }

    cout << excluded_count << endl;

    return 0;
}