#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

const int MOD = 1e9 + 7;
const int MAX = 6005;

long long fact[MAX], invFact[MAX], pow2[MAX], invPow2[MAX];

// 분할 정복을 이용한 거듭제곱
long long power(long long base, long long exp) {
    long long res = 1;
    base %= MOD;
    while (exp > 0) {
        if (exp % 2 == 1) res = (res * base) % MOD;
        base = (base * base) % MOD;
        exp /= 2;
    }
    return res;
}

// 모듈로 역원 계산 (페르마의 소정리)
long long modInverse(long long n) {
    return power(n, MOD - 2);
}

// 팩토리얼 및 역원 사전 계산
void precompute() {
    fact[0] = 1;
    invFact[0] = 1;
    pow2[0] = 1;
    invPow2[0] = 1;
    
    long long inv2 = modInverse(2);
    for (int i = 1; i < MAX; i++) {
        fact[i] = (fact[i - 1] * i) % MOD;
        invFact[i] = modInverse(fact[i]);
        pow2[i] = (pow2[i - 1] * 2) % MOD;
        invPow2[i] = (invPow2[i - 1] * inv2) % MOD;
    }
}

// 조합 (nCr) 계산 함수
long long nCr(int n, int r) {
    if (r < 0 || r > n) return 0;
    return fact[n] * invFact[r] % MOD * invFact[n - r] % MOD;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    precompute();
    
    int n, m;
    if (!(cin >> n >> m)) return 0;
    
    long long ans = 0;
    
    int min_S = max(n, m);
    int max_S = min(2 * n, 2 * m);
    
    // 가능한 1의 총 개수 S에 대해 반복
    for (int s = min_S; s <= max_S; s++) {
        int r1 = 2 * n - s;
        int r2 = s - n;
        int c1 = 2 * m - s;
        int c2 = s - m;
        
        // 차수가 1, 2인 행과 열을 배치하는 경우의 수
        long long ways_choose = nCr(n, r1) * nCr(m, c1) % MOD;
        
        long long inner_sum = 0;
        int limit = min(r2, c2);
        
        // 포함-배제의 원리를 이용해 이중 간선(Double Edge)이 없는 소켓 매칭 계산
        for (int k = 0; k <= limit; k++) {
            long long term = nCr(r2, k) * nCr(c2, k) % MOD;
            term = term * fact[k] % MOD;
            term = term * pow2[k] % MOD;
            term = term * fact[s - 2 * k] % MOD;
            
            if (k % 2 == 1) {
                inner_sum = (inner_sum - term + MOD) % MOD;
            } else {
                inner_sum = (inner_sum + term) % MOD;
            }
        }
        
        // 같은 정점의 두 소켓 중복 제거 (2^(r2 + c2) 로 나누기)
        long long ways_arrange = inner_sum * invPow2[r2 + c2] % MOD;
        
        // 정답 누적
        ans = (ans + ways_choose * ways_arrange) % MOD;
    }
    
    cout << ans << "\n";
    return 0;
}