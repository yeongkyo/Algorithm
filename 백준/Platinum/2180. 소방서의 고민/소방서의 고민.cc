#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

struct Fire {
    long long a, b;
    bool operator<(const Fire& other) const {
        if (a == 0 && other.a == 0) return b < other.b;
        if (a == 0) return false;
        if (other.a == 0) return true;
        if (b * other.a != other.b * a) {
            return b * other.a < other.b * a;
        }
        if (a != other.a) return a < other.a;
        return b < other.b;
    }
};

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    if (!(cin >> n)) return 0;

    vector<Fire> fires(n);
    for (int i = 0; i < n; ++i) {
        cin >> fires[i].a >> fires[i].b;
    }

    sort(fires.begin(), fires.end());

    long long T = 0;
    for (int i = 0; i < n; ++i) {
        T = ((fires[i].a + 1) * T + fires[i].b) % 40000;
    }

    cout << T << "\n";

    return 0;
}