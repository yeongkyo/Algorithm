#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

struct Result {
    bool valid;
    string A, B, C;
    bool operator<(const Result& other) const {
        if (!valid) return true;
        if (!other.valid) return false;
        if (C != other.C) return C < other.C;
        if (A != other.A) return A < other.A;
        return false;
    }
};

string sA, sB, sC, rA, rB, rC;
int lenA, lenB, lenC, N;
bool visited[55][2];
Result memo[55][2];

vector<int> get_valid(char c, int idx, int L) {
    if (idx >= L) return {0};
    vector<int> res;
    if (c == '?') {
        for (int d = 0; d <= 9; ++d) {
            if (d == 0 && L > 1 && idx == L - 1) continue;
            res.push_back(d);
        }
    } else {
        int d = c - '0';
        if (d == 0 && L > 1 && idx == L - 1) return {}; 
        res.push_back(d);
    }
    return res;
}

Result solve(int idx, int carry) {
    if (idx == N) {
        if (carry == 0) return {true, "", "", ""};
        else return {false, "", "", ""};
    }
    if (visited[idx][carry]) return memo[idx][carry];
    
    Result best = {false, "", "", ""};
    
    vector<int> cA = get_valid(idx < lenA ? rA[idx] : ' ', idx, lenA);
    vector<int> cB = get_valid(idx < lenB ? rB[idx] : ' ', idx, lenB);
    vector<int> cC = get_valid(idx < lenC ? rC[idx] : ' ', idx, lenC);
    
    for (int dA : cA) {
        for (int dB : cB) {
            int sum = dA + dB + carry;
            int dC = sum % 10;
            int n_carry = sum / 10;
            
            bool valid_C = false;
            for (int x : cC) {
                if (x == dC) {
                    valid_C = true;
                    break;
                }
            }
            
            if (valid_C) {
                Result nxt = solve(idx + 1, n_carry);
                if (nxt.valid) {
                    Result cur;
                    cur.valid = true;
                    cur.A = nxt.A; if (idx < lenA) cur.A += (char)(dA + '0');
                    cur.B = nxt.B; if (idx < lenB) cur.B += (char)(dB + '0');
                    cur.C = nxt.C; if (idx < lenC) cur.C += (char)(dC + '0');
                    if (best < cur) {
                        best = cur;
                    }
                }
            }
        }
    }
    
    visited[idx][carry] = true;
    return memo[idx][carry] = best;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    string s;
    if (cin >> s) {
        int plus = s.find('+');
        int eq = s.find('=');
        
        sA = s.substr(0, plus);
        sB = s.substr(plus + 1, eq - plus - 1);
        sC = s.substr(eq + 1);
        
        rA = sA; reverse(rA.begin(), rA.end());
        rB = sB; reverse(rB.begin(), rB.end());
        rC = sC; reverse(rC.begin(), rC.end());
        
        lenA = rA.length();
        lenB = rB.length();
        lenC = rC.length();
        
        N = max({lenA, lenB, lenC});
        
        Result res = solve(0, 0);
        
        if (res.valid) cout << res.A << "+" << res.B << "=" << res.C << "\n";
        else cout << "-1\n";
    }
    return 0;
}