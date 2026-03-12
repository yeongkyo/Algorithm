#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

struct Route {
    long long start;
    long long end;
    int id;
    bool operator<(const Route& other) const {
        if (start != other.start)
            return start < other.start;
        return end > other.end;
    }
};

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    long long N;
    int M;
    if (!(cin >> N >> M)) return 0;

    vector<Route> routes;
    routes.reserve(2 * M);

    for (int i = 1; i <= M; ++i) {
        long long a, b;
        cin >> a >> b;
        if (a < b) {
            routes.push_back({a, b, i});
            routes.push_back({a + N, b + N, i});
        } else {
            routes.push_back({a, b + N, i});
        }
    }

    sort(routes.begin(), routes.end());

    vector<bool> is_contained(M + 1, false);
    long long max_end = -1;

    for (const auto& route : routes) {
        if (route.end <= max_end) {
            is_contained[route.id] = true;
        } else {
            max_end = route.end;
        }
    }

    for (int i = 1; i <= M; ++i) {
        if (!is_contained[i]) {
            cout << i << " ";
        }
    }
    cout << "\n";

    return 0;
}