#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    cin >> n;

    long long a, b;
    cin >> a >> b;

    vector<pair<long long, long long>> points(n);
    for (int i = 0; i < n; i++) {
        cin >> points[i].first >> points[i].second;
    }

    sort(points.begin(), points.end());

    int count = 0;
    for (int i = 0; i < n; i++) {
        long long x = points[i].first;
        long long y = points[i].second;

        if (binary_search(points.begin(), points.end(), make_pair(x + a, y)) &&
            binary_search(points.begin(), points.end(), make_pair(x, y + b)) &&
            binary_search(points.begin(), points.end(), make_pair(x + a, y + b))) {
            count++;
        }
    }

    cout << count << "\n";

    return 0;
}