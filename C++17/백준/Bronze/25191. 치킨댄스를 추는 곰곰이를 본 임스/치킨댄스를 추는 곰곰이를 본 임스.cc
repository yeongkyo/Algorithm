#include <iostream>
#include <algorithm>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int N, A, B;
    cin >> N;
    cin >> A >> B;

    int max_with_drinks = (A / 2) + B;

    cout << min(N, max_with_drinks) << "\n";

    return 0;
}