#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int TC;
    cin >> TC;
    while (TC--) {
        int N;
        string s;
        cin >> N >> s;

        string st;
        st.reserve(N);

        for (char c : s) {
            st.push_back(c);
            if (st.size() >= 3) {
                int m = (int)st.size();
                if (st[m-3] == 'f' && st[m-2] == 'o' && st[m-1] == 'x') {
                    st.resize(m - 3);
                }
            }
        }

        cout << st.size() << '\n';
    }
    return 0;
}