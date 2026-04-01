#include <iostream>
#include <algorithm>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    long long N, M;
    cin >> N >> M;

    long long s_sets = N / 2;
    long long a_sets = M / 2;

    cout << min(s_sets, a_sets) << "\n";

    return 0;
}