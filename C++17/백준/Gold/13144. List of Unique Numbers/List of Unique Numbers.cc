#include <iostream>
#include <vector>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    if (!(cin >> n)) return 0;

    vector<int> arr(n);
    for (int i = 0; i < n; i++) {
        cin >> arr[i];
    }

    vector<bool> visited(100001, false);
    long long count = 0;
    int end = 0;

    for (int start = 0; start < n; start++) {
        while (end < n && !visited[arr[end]]) {
            visited[arr[end]] = true;
            end++;
        }
        
        count += (end - start);
        visited[arr[start]] = false;
    }

    cout << count << "\n";

    return 0;
}