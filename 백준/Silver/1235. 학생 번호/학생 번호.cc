#include <iostream>
#include <vector>
#include <string>
#include <set>

using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    cin >> n;
    vector<string> v(n);
    for (int i = 0; i < n; i++) {
        cin >> v[i];
    }

    int len = v[0].length();
    for (int k = 1; k <= len; k++) {
        set<string> s;
        bool is_unique = true;
        for (int i = 0; i < n; i++) {
            string suffix = v[i].substr(len - k);
            if (s.count(suffix)) {
                is_unique = false;
                break;
            }
            s.insert(suffix);
        }
        if (is_unique) {
            cout << k << "\n";
            break;
        }
    }

    return 0;
}