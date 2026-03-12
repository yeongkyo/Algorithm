#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
using namespace std;

typedef unsigned long long ull;

// forwardDP_bp: 
// s[a1,a2)와 t[b1,b2)에 대해 비트병렬 LCS 계산을 하여
// F[j] = LCS( s-substring, t[b1, b1+j] ) (j=0..(b2-b1)) 를 반환한다.
vector<int> forwardDP_bp(const string &s, int a1, int a2,
                         const string &t, int b1, int b2) {
    int n = a2 - a1;   // s 구간 길이
    int m = b2 - b1;   // t 구간 길이
    // s 구간을 64비트 블록으로 나눈 개수
    int blocks = (n + 63) / 64;
    // 미리 s 구간에 대한 각 문자('A'-'Z')의 bitmask를 계산
    vector< vector<ull> > mask(26, vector<ull>(blocks, 0ULL));
    for (int i = 0; i < n; i++) {
        int c = s[a1 + i] - 'A';
        int bi = i / 64;
        int pos = i % 64;
        mask[c][bi] |= (1ULL << pos);
    }
    
    // dp: s 구간에 대해 현재 매칭 상태를 bitset으로 저장 (초기 모두 0)
    vector<ull> dp(blocks, 0ULL);
    // F[j] : t[b1, b1+j] 까지 처리한 후 LCS 길이 (bit count)
    vector<int> F(m + 1, 0);
    F[0] = 0;
    
    for (int j = 0; j < m; j++) {
        int c = t[b1 + j] - 'A';
        // X = dp OR mask[c]
        vector<ull> X(blocks, 0ULL);
        for (int i = 0; i < blocks; i++) {
            X[i] = dp[i] | mask[c][i];
        }
        // Y = (dp << 1) OR 1, (블록간 carry 고려)
        vector<ull> Y(blocks, 0ULL);
        ull carry = 0;
        for (int i = 0; i < blocks; i++) {
            ull shifted = (dp[i] << 1) | carry;
            Y[i] = shifted;
            carry = (dp[i] >> 63) & 1ULL;
        }
        Y[0] |= 1ULL;  // 전체 bitset의 최하위 비트 1
        if(n % 64 != 0)
            Y[blocks-1] &= ((1ULL << (n % 64)) - 1);
        
        // sub = X - Y (여러 블록에 대해 borrow 고려)
        vector<ull> sub(blocks, 0ULL);
        ull bor = 0;
        for (int i = 0; i < blocks; i++) {
            ull tmp = X[i] - Y[i] - bor;
            if (X[i] < Y[i] + bor) bor = 1;
            else bor = 0;
            sub[i] = tmp;
        }
        // dp = X & ~(sub)
        for (int i = 0; i < blocks; i++) {
            dp[i] = X[i] & ~(sub[i]);
        }
        // F[j+1] = dp의 전체 1의 개수 = LCS 길이까지의 값
        int cnt = 0;
        for (int i = 0; i < blocks; i++) {
            cnt += __builtin_popcountll(dp[i]);
        }
        F[j+1] = cnt;
    }
    return F;
}

// backwardDP_bp: s[a1,a2), t[b1,b2)에 대해 후방(dp역방향) 값을 계산한다.
// s와 t의 해당 구간을 뒤집어 forwardDP_bp를 호출한 후, 결과 배열을 뒤집어서 반환한다.
vector<int> backwardDP_bp(const string &s, int a1, int a2,
                          const string &t, int b1, int b2) {
    int n = a2 - a1, m = b2 - b1;
    string s_rev, t_rev;
    s_rev.resize(n); t_rev.resize(m);
    for (int i = 0; i < n; i++)
        s_rev[i] = s[a2 - 1 - i];
    for (int i = 0; i < m; i++)
        t_rev[i] = t[b2 - 1 - i];
    
    vector<int> R = forwardDP_bp(s_rev, 0, n, t_rev, 0, m);
    // R: R[0..m] 에 대해, R[m - j]가 우리가 필요한 후방값 G[j]
    vector<int> G(m + 1, 0);
    for (int j = 0; j <= m; j++) {
        G[j] = R[m - j];
    }
    return G;
}

// hirschberg_bp: 비트병렬 forward/backward DP를 이용한 Hirschberg 재귀적 LCS 복원 함수
string hirschberg_bp(const string &s, int a1, int a2,
                     const string &t, int b1, int b2) {
    if (a1 >= a2) return "";
    // a 구간 길이가 1이면, s[a1]가 t 구간 내에 존재하면 그 문자 반환
    if (a2 - a1 == 1) {
        for (int j = b1; j < b2; j++) {
            if (s[a1] == t[j])
                return string(1, s[a1]);
        }
        return "";
    }
    
    int aMid = (a1 + a2) / 2;
    vector<int> L1 = forwardDP_bp(s, a1, aMid, t, b1, b2);
    vector<int> L2 = backwardDP_bp(s, aMid, a2, t, b1, b2);
    int m = b2 - b1;
    int partition = 0, best = -1;
    // 적절한 B분할 지점을 찾는다.
    for (int j = 0; j <= m; j++) {
        int sum = L1[j] + L2[j];
        if (sum > best) {
            best = sum;
            partition = j;
        }
    }
    string left = hirschberg_bp(s, a1, aMid, t, b1, b1 + partition);
    string right = hirschberg_bp(s, aMid, a2, t, b1 + partition, b2);
    return left + right;
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    
    string s1, s2;
    cin >> s1 >> s2;
    
    string lcs = hirschberg_bp(s1, 0, s1.size(), s2, 0, s2.size());
    cout << lcs.size() << "\n" << lcs << "\n";
    return 0;
}
