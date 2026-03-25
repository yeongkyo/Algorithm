#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

typedef long long ll;

int main() {
    // 입출력 속도 최적화
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int N, K;
    if (!(cin >> N >> K)) return 0;

    vector<ll> A(N + 1);
    vector<ll> P(N + 1, 0); // 누적 합 (Prefix Sum)
    for (int i = 1; i <= N; i++) {
        cin >> A[i];
        P[i] = P[i - 1] + A[i];
    }

    vector<ll> L(N + 1, 0);
    vector<ll> M(N + 1, 0);

    L[0] = 0;
    M[0] = 0; // L[0] - P[0]

    for (int i = 1; i < N; i++) {
        // 1. i번째 아이템을 줍지 않는 경우
        L[i] = L[i - 1];

        // 2. i번째 아이템이 1번과 연결된 접두사의 일부인 경우
        L[i] = max(L[i], P[i]);

        // 3. i번째 아이템이 중간 블록(길이 >= K)의 끝인 경우
        // j+1부터 i까지의 블록을 선택한다고 하면 (i - j >= K)
        // L[i] = max(L[j] + P[i] - P[j]) = P[i] + max(L[j] - P[j])
        if (i >= K) {
            L[i] = max(L[i], P[i] + M[i - K]);
        }

        // 다음 단계 계산을 위해 M[i] 갱신
        M[i] = max(M[i - 1], L[i] - P[i]);
    }
    ll ans = L[N - 1]; // N을 줍지 않는 경우
    ans = max(ans, P[N]); // 전체를 다 줍는 경우
    ans = max(ans, P[N] + M[N - 1]); // 임의의 j까지 최적으로 줍고, [j+1, N]을 접미사로 줍는 경우

    cout << ans << "\n";

    return 0;
}