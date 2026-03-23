#include <iostream>
#include <stack>
#include <vector>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    cin >> n;

    stack<int> s;
    int count = 0;

    for (int i = 0; i < n; ++i) {
        int x, y;
        cin >> x >> y;

        while (!s.empty() && s.top() > y) {
            count++;
            s.pop();
        }

        if (!s.empty() && s.top() == y) {
            continue;
        }

        if (y > 0) {
            s.push(y);
        }
    }

    while (!s.empty()) {
        count++;
        s.pop();
    }

    cout << count << "\n";

    return 0;
}