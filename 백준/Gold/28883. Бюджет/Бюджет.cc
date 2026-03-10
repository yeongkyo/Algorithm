#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

using namespace std;

// 잔고 범위가 -100 ~ 100이므로 안전하게 100을 더해 인덱스로 사용
const int OFFSET = 100;
bool dp[10001][205];
int choice[10001][205]; // 1: 수입, 0: 지출

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    int n, a, b;
    cin >> n >> a >> b;

    vector<int> f(n + 1);
    for (int i = 1; i <= n; i++) {
        cin >> f[i];
    }

    // 초기 상태: 0번째 연산 후 잔고는 0
    dp[0][0 + OFFSET] = true;

    for (int i = 1; i <= n; i++) {
        for (int j = a; j <= b; j++) {
            // 1. 수입인 경우 (이전 잔고 + f[i] = 현재 잔고 j)
            int prev_income = j - f[i];
            // 0번째 잔고는 a, b 범위 제한이 없으므로 별도 처리하거나 범위를 넉넉히 잡음
            if (i == 1) {
                if (prev_income == 0) {
                    dp[i][j + OFFSET] = true;
                    choice[i][j + OFFSET] = 1;
                }
            } else if (prev_income >= a && prev_income <= b) {
                if (dp[i - 1][prev_income + OFFSET]) {
                    dp[i][j + OFFSET] = true;
                    choice[i][j + OFFSET] = 1;
                }
            }

            // 2. 지출인 경우 (이전 잔고 - f[i] = 현재 잔고 j)
            if (dp[i][j + OFFSET]) continue; // 이미 수입으로 가능하면 넘어감 (아무거나 출력해도 되므로)

            int prev_expense = j + f[i];
            if (i == 1) {
                if (prev_expense == 0) {
                    dp[i][j + OFFSET] = true;
                    choice[i][j + OFFSET] = 0;
                }
            } else if (prev_expense >= a && prev_expense <= b) {
                if (dp[i - 1][prev_expense + OFFSET]) {
                    dp[i][j + OFFSET] = true;
                    choice[i][j + OFFSET] = 0;
                }
            }
        }
    }

    // 결과 확인 및 역추적
    int last_balance = -1000;
    for (int j = a; j <= b; j++) {
        if (dp[n][j + OFFSET]) {
            last_balance = j;
            break;
        }
    }

    if (last_balance == -1000) {
        cout << "Impossible" << endl;
    } else {
        string result = "";
        int current_j = last_balance;
        for (int i = n; i >= 1; i--) {
            int c = choice[i][current_j + OFFSET];
            result += to_string(c);
            if (c == 1) current_j -= f[i]; // 수입이었으면 이전엔 더 작았음
            else current_j += f[i];       // 지출이었으면 이전엔 더 컸음
        }
        reverse(result.begin(), result.end());
        cout << result << endl;
    }

    return 0;
}