#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int N;
    char game;
    cin >> N >> game;

    unordered_set<string> s;
    s.reserve(N * 2);

    for (int i = 0; i < N; i++) {
        string name;
        cin >> name;
        s.insert(name);
    }

    int need = (game == 'Y' ? 1 : (game == 'F' ? 2 : 3));
    cout << (int)s.size() / need << "\n";
    return 0;
}