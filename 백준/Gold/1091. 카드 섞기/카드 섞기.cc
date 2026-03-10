#include <iostream>
#include <vector>

using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    cin >> n;

    vector<int> p(n);
    for (int i = 0; i < n; i++) cin >> p[i];

    vector<int> s(n);
    for (int i = 0; i < n; i++) cin >> s[i];

    vector<int> cards(n);
    vector<int> initial_cards(n);
    for (int i = 0; i < n; i++) {
        cards[i] = i;
        initial_cards[i] = i;
    }

    int count = 0;
    while (true) {
        bool match = true;
        for (int i = 0; i < n; i++) {
            if (p[cards[i]] != i % 3) {
                match = false;
                break;
            }
        }

        if (match) {
            cout << count << "\n";
            return 0;
        }

        vector<int> next_cards(n);
        for (int i = 0; i < n; i++) {
            next_cards[s[i]] = cards[i];
        }
        cards = next_cards;
        count++;

        bool is_initial = true;
        for (int i = 0; i < n; i++) {
            if (cards[i] != initial_cards[i]) {
                is_initial = false;
                break;
            }
        }
        if (is_initial) break;
    }

    cout << -1 << "\n";

    return 0;
}