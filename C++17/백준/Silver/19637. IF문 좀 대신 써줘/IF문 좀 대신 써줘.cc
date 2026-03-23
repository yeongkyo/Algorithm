#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

struct Title {
    string name;
    int limit;
};

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n, m;
    cin >> n >> m;

    vector<Title> titles;
    for (int i = 0; i < n; i++) {
        string name;
        int limit;
        cin >> name >> limit;
        if (!titles.empty() && titles.back().limit == limit) continue;
        titles.push_back({name, limit});
    }

    for (int i = 0; i < m; i++) {
        int power;
        cin >> power;

        int low = 0, high = (int)titles.size() - 1;
        int ans_idx = 0;

        while (low <= high) {
            int mid = (low + high) / 2;
            if (titles[mid].limit >= power) {
                ans_idx = mid;
                high = mid - 1;
            } else {
                low = mid + 1;
            }
        }
        cout << titles[ans_idx].name << "\n";
    }

    return 0;
}