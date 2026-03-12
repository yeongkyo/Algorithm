#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

typedef __int128_t int128;

struct Cand {
    int128 len;
    int128 pos;
};

int128 B1 = 1000000;

int128 get_Offset(int128 N, int128 L1, int128 L2) {
    if (N <= 0) return 0;
    return (N - 1) * B1 * L1 + (N - 1) * N / 2 * L2;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    string S1, S2;
    long long c_in;
    if (!(cin >> S1 >> S2 >> c_in)) return 0;
    int128 C = c_in;

    int128 L1 = S1.length();
    int128 L2 = S2.length();

    int128 P1 = 0, Suf1 = 0;
    while (P1 < L1 && S1[P1] == '0') P1++;
    if (P1 == L1) {
        Suf1 = L1;
    } else {
        while (Suf1 < L1 && S1[L1 - 1 - Suf1] == '0') Suf1++;
    }
    bool S1_all = (P1 == L1);

    int128 P2 = 0, Suf2 = 0;
    while (P2 < L2 && S2[P2] == '0') P2++;
    if (P2 == L2) {
        Suf2 = L2;
    } else {
        while (Suf2 < L2 && S2[L2 - 1 - Suf2] == '0') Suf2++;
    }
    bool S2_all = (P2 == L2);

    vector<Cand> cands;

    auto add_z_runs = [&](const string& s, int128 offset) {
        int128 cur = 0;
        for (int i = 0; i < s.length(); ++i) {
            if (s[i] == '0') {
                cur++;
            } else {
                if (cur > 0) {
                    cands.push_back({cur, offset + i - cur});
                    cur = 0;
                }
            }
        }
        if (cur > 0) {
            cands.push_back({cur, offset + s.length() - cur});
        }
    };

    if (!S1_all && !S2_all) {
        add_z_runs(S1, 0);
        cands.push_back({Suf1 + P1, L1 - Suf1});
        cands.push_back({Suf1 + P2, B1 * L1 - Suf1});
        add_z_runs(S2, B1 * L1);
        cands.push_back({Suf2 + P2, get_Offset(2, L1, L2) + B1 * L1 + L2 - Suf2});
        cands.push_back({Suf2 + P1, get_Offset(2, L1, L2) - Suf2});
    }
    else if (S1_all && !S2_all) {
        cands.push_back({B1 * L1 + P2, 0});
        add_z_runs(S2, B1 * L1);
        cands.push_back({Suf2 + P2, get_Offset(2, L1, L2) + B1 * L1 + L2 - Suf2});
        cands.push_back({Suf2 + B1 * L1 + P2, get_Offset(2, L1, L2) - Suf2});
    }
    else if (!S1_all && S2_all) {
        add_z_runs(S1, 0);
        cands.push_back({Suf1 + P1, L1 - Suf1});
        
        int128 req_N = 1;
        if (C > Suf1 + P1) {
            req_N = (C - Suf1 - P1 + L2 - 1) / L2;
            if (req_N < 1) req_N = 1;
        }
        cands.push_back({Suf1 + req_N * L2 + P1, get_Offset(req_N, L1, L2) + B1 * L1 - Suf1});
    }
    else if (S1_all && S2_all) {
        int128 inf = 20000000000000000LL; 
        cands.push_back({inf, 0});
    }

    int128 MAX_POS = 10000000000000000LL;
    int128 min_pos = -1;

    for (auto& cand : cands) {
        if (cand.len >= C) {
            if (cand.pos + C <= MAX_POS) {
                if (min_pos == -1 || cand.pos < min_pos) {
                    min_pos = cand.pos;
                }
            }
        }
    }

    if (min_pos == -1) {
        cout << -1 << "\n";
    } else {
        long long ans = (long long)min_pos;
        cout << ans << "\n";
    }

    return 0;
}