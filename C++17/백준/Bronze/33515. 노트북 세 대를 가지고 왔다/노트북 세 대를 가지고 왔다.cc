#include <iostream>
#include <algorithm>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int t1, t2;
    cin >> t1 >> t2;

    cout << min(t1, t2) << "\n";

    return 0;
}