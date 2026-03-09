#include <iostream>
#include <vector>

using namespace std;

void solve() {
    int n;
    cin >> n;
    while (n--) {
        int a_cnt[5] = {0}, b_cnt[5] = {0};
        int num, shape;

        // 어린이 A 정보 입력
        cin >> num;
        for (int i = 0; i < num; ++i) {
            cin >> shape;
            a_cnt[shape]++;
        }

        // 어린이 B 정보 입력
        cin >> num;
        for (int i = 0; i < num; ++i) {
            cin >> shape;
            b_cnt[shape]++;
        }

        // 우선순위가 높은 4(별)부터 1(세모)까지 비교
        bool settled = false;
        for (int i = 4; i >= 1; --i) {
            if (a_cnt[i] > b_cnt[i]) {
                cout << "A\n";
                settled = true;
                break;
            } else if (a_cnt[i] < b_cnt[i]) {
                cout << "B\n";
                settled = true;
                break;
            }
        }

        if (!settled) cout << "D\n";
    }
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    solve();
    return 0;
}