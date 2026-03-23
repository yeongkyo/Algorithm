#include <iostream>
#include <string>
#include <unordered_set>
#include <sstream>

using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    int n, m;
    cin >> n >> m;

    unordered_set<string> s;
    for (int i = 0; i < n; i++) {
        string word;
        cin >> word;
        s.insert(word);
    }

    for (int i = 0; i < m; i++) {
        string line;
        cin >> line;
        stringstream ss(line);
        string keyword;

        while (getline(ss, keyword, ',')) {
            s.erase(keyword);
        }
        cout << s.size() << "\n";
    }

    return 0;
}