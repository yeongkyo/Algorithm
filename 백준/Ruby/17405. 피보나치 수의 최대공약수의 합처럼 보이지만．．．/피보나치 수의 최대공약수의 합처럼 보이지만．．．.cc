#include <iostream>
#include <vector>

using namespace std;

const long long MOD = 1000000007;
const int M = 2000000;

int min_prime[M + 1];
vector<int> primes;
int phi_val[M + 1];
long long Phi[M + 1];

// 시간 초과 해결의 핵심: x^k, F_x, P_H(x)를 선형 체와 함께 한 번에 캐싱
long long pk[M + 1]; 
long long F[M + 1];
long long PH[M + 1];

long long n, k;

// 거듭제곱 계산 (분할 정복)
long long power(long long base, long long exp) {
    long long res = 1;
    base %= MOD;
    while (exp > 0) {
        if (exp % 2 == 1) res = res * base % MOD;
        base = base * base % MOD;
        exp /= 2;
    }
    return res;
}

// 확장된 선형 체 알고리즘
void sieve() {
    phi_val[1] = 1;
    pk[1] = 1; // 1^k = 1
    for (int i = 2; i <= M; ++i) {
        if (min_prime[i] == 0) {
            min_prime[i] = i;
            primes.push_back(i);
            phi_val[i] = i - 1;
            pk[i] = power(i, k); // 소수에 대해서만 거듭제곱 연산
        }
        for (int p : primes) {
            if (p > min_prime[i] || i * p > M) break;
            min_prime[i * p] = p;
            pk[i * p] = pk[i] * pk[p] % MOD; // 곱셈적 성질 활용
            if (i % p == 0) {
                phi_val[i * p] = phi_val[i] * p;
            } else {
                phi_val[i * p] = phi_val[i] * (p - 1);
            }
        }
    }
    
    // 배열 초기값 세팅 및 누적합/피보나치 점화식 계산
    Phi[0] = 0;
    PH[0] = 0;
    Phi[1] = 1;
    F[1] = 1;
    if (M >= 2) F[2] = 1;
    PH[1] = pk[1] * F[1] % MOD;
    
    for (int i = 2; i <= M; ++i) {
        if (i > 2) F[i] = (F[i - 1] + F[i - 2]) % MOD;
        PH[i] = (PH[i - 1] + pk[i] * F[i]) % MOD;
        Phi[i] = (Phi[i - 1] + phi_val[i]) % MOD;
    }
}

long long memo_arr[2005];
bool visited[2005];

// 두 체 (Du Sieve)를 이용한 대규모 Phi 계산
long long get_Phi(long long v) {
    if (v <= M) return Phi[v];
    int idx = n / v; 
    if (visited[idx]) return memo_arr[idx];

    long long res = v % MOD * ((v + 1) % MOD) % MOD * 500000004 % MOD;
    for (long long l = 2, r; l <= v; l = r + 1) {
        r = v / (v / l);
        long long term = (r - l + 1) % MOD;
        res = (res - term * get_Phi(v / l)) % MOD;
    }
    res = (res + MOD) % MOD;
    visited[idx] = true;
    return memo_arr[idx] = res;
}

// 조합(nCr) 전처리 상수
long long fact[4005], invFact[4005];

long long modInverse(long long x) {
    return power(x, MOD - 2);
}

void precompute_nCr() {
    fact[0] = 1;
    invFact[0] = 1;
    for (int i = 1; i <= 4000; i++) {
        fact[i] = (fact[i - 1] * i) % MOD;
    }
    invFact[4000] = modInverse(fact[4000]);
    for (int i = 3999; i >= 1; i--) {
        invFact[i] = (invFact[i + 1] * (i + 1)) % MOD;
    }
}

long long nCr(int n_val, int r_val) {
    if (r_val < 0 || r_val > n_val) return 0;
    return fact[n_val] * invFact[r_val] % MOD * invFact[n_val - r_val] % MOD;
}

long long a_poly[4005];

// 호너의 법칙(Horner's Method)을 이용한 다항식 A(x) 빠른 계산
long long eval_A(long long x) {
    x %= MOD;
    long long res = 0;
    for (int i = k; i >= 0; --i) {
        res = (res * x + a_poly[i]) % MOD;
    }
    return res;
}

// 행렬 거듭제곱을 응용한 O(log x) 피보나치 계산
pair<long long, long long> fib(long long x) {
    if (x == 0) return {0, 1};
    auto p = fib(x >> 1);
    long long c = p.first * (2 * p.second - p.first + MOD) % MOD;
    long long d = (p.first * p.first + p.second * p.second) % MOD;
    if (x & 1) return {d, (c + d) % MOD};
    else return {c, d};
}

// 부분합 P_H(x) 구하기 (가장 중요한 캐싱 로직 적용 부분)
long long get_P_H(long long x) {
    if (x <= M) return PH[x]; // M 이하라면 전처리된 값을 즉시 반환 O(1)
    
    // M을 초과하는 매우 큰 경우에만 다항식 직접 연산 수행
    long long Ax = eval_A(x);
    long long Ax1 = eval_A(x + 1);
    auto f = fib(x);
    long long Fx = f.first;
    long long Fx1 = f.second;
    long long A0 = a_poly[0];
    
    long long res = (Ax * Fx1 % MOD + Ax1 * Fx % MOD - A0) % MOD;
    return (res + MOD) % MOD;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    if (!(cin >> n >> k)) return 0;

    sieve();
    precompute_nCr();

    // 다항식 A(x)의 계수들을 전처리
    a_poly[k] = 1;
    for (int m = k - 1; m >= 0; --m) {
        long long sum = 0;
        for (int j = m + 1; j <= k; j += 2) {
            sum = (sum + 2 * nCr(j, m) % MOD * a_poly[j]) % MOD;
        }
        a_poly[m] = (MOD - sum) % MOD;
    }

    long long ans = 0;
    
    // 블록 단위 합산
    for (long long l = 1, r; l <= n; l = r + 1) {
        long long v = n / l;
        r = n / v;
        
        long long f_v = (2 * get_Phi(v) - 1) % MOD;
        f_v = (f_v + MOD) % MOD;
        
        long long PH_r = get_P_H(r);
        long long PH_l_1 = get_P_H(l - 1);
        long long term = (PH_r - PH_l_1 + MOD) % MOD;
        
        ans = (ans + f_v * term) % MOD;
    }

    cout << ans << "\n";
    return 0;
}