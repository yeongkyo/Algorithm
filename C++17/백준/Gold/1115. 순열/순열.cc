#include <iostream>
#include <vector>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int N;
    if (!(cin >> N)) return 0;

    vector<int> P(N);
    for (int i = 0; i < N; ++i) {
        cin >> P[i];
    }

    vector<bool> visited(N, false);
    int cycle_count = 0;

    for (int i = 0; i < N; ++i) {
        if (!visited[i]) {
            cycle_count++;
            int curr = i;
            while (!visited[curr]) {
                visited[curr] = true;
                curr = P[curr];
            }
        }
    }

    if (cycle_count == 1) {
        cout << 0 << "\n";
    } else {
        cout << cycle_count << "\n";
    }

    return 0;
}