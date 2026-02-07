#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

// DP 상태를 저장할 구조체
struct Node {
    long long not_picked; // 해당 노드를 선택하지 않았을 때의 최대 신뢰도
    long long picked;     // 해당 노드를 선택했을 때의 최대 신뢰도
};

int main() {
    // 입출력 속도 최적화
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    if (!(cin >> n)) return 0;

    // 신뢰도 입력
    vector<long long> w(n);
    for (int i = 0; i < n; ++i) {
        cin >> w[i];
    }

    // 단계별 정보 입력
    // host[i]와 protocol[i]는 i번 노드가 추가될 때의 정보를 저장
    vector<int> host(n), protocol(n);
    for (int i = 1; i < n; ++i) {
        cin >> host[i] >> protocol[i];
    }

    // DP 테이블 초기화
    vector<Node> dp(n);
    for (int i = 0; i < n; ++i) {
        dp[i].not_picked = 0;
        dp[i].picked = w[i];
    }

    // 역순으로 처리 (가장 나중에 추가된 노드부터 병합)
    for (int i = n - 1; i >= 1; --i) {
        int h = host[i];
        int p = protocol[i];

        if (p == 0) {
            // Protocol 0: IAmYourFriend (Leaf node logic)
            // h가 선택되면 i는 선택 불가
            dp[h].picked += dp[i].not_picked;
            // h가 선택되지 않으면 i는 선택되어도 되고 안 되어도 됨
            dp[h].not_picked += max(dp[i].picked, dp[i].not_picked);
        } 
        else if (p == 1) {
            // Protocol 1: MyFriendsAreYourFriends (False Twins)
            // h가 선택되면 이웃들이 차단됨. i는 이웃들과만 연결되어 있으므로 i는 자유로움.
            // i와 h 사이에는 간선이 없으므로 h가 선택되어도 i 선택 가능.
            dp[h].picked += max(dp[i].picked, dp[i].not_picked);
            // h가 선택되지 않은 상태(이웃이 선택될 수도 있는 상태)라면 i도 이웃에 의해 차단될 수 있으므로 선택 불가 상태 합산
            dp[h].not_picked += dp[i].not_picked;
        } 
        else if (p == 2) {
            // Protocol 2: WeAreYourFriends (True Twins)
            // h와 i는 서로 연결되어 있고 이웃도 공유함.
            // 'Picked' 상태: 둘 중 하나만 선택 가능
            // 1) h 선택, i 선택 안 함: dp[h].picked + dp[i].not_picked
            // 2) i 선택, h 선택 안 함: dp[h].not_picked + dp[i].picked 
            // (i를 선택하는 것은 h를 선택하는 것과 외부적으로 동일한 효과를 가짐)
            dp[h].picked = max(dp[h].picked + dp[i].not_picked, dp[h].not_picked + dp[i].picked);
            
            // 'Not Picked' 상태: 둘 다 선택 안 함
            dp[h].not_picked += dp[i].not_picked;
        }
    }

    // 최종 결과: 노드 0(루트)을 선택했을 때와 안 했을 때 중 최댓값
    cout << max(dp[0].picked, dp[0].not_picked) << endl;

    return 0;
}