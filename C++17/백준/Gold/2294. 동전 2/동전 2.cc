#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

int main() {
    // 입출력 속도 향상
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n, k;
    cin >> n >> k;

    vector<int> coins(n);
    for (int i = 0; i < n; i++) {
        cin >> coins[i];
    }

    // dp[i]를 k+1로 초기화 (최대 k개까지 사용 가능하므로 k+1은 INF의 의미)
    vector<int> dp(k + 1, 10001);
    dp[0] = 0;

    for (int i = 0; i < n; i++) {
        int current_coin = coins[i];
        // 현재 동전의 가치부터 목표 금액 k까지 업데이트
        for (int j = current_coin; j <= k; j++) {
            // j원을 만드는 데 현재 동전을 사용하는 게 이득인지 확인
            if (dp[j - current_coin] != 10001) {
                dp[j] = min(dp[j], dp[j - current_coin] + 1);
            }
        }
    }

    // 결과 출력
    if (dp[k] == 10001) {
        cout << -1 << "\n"; // 불가능한 경우
    } else {
        cout << dp[k] << "\n";
    }

    return 0;
}