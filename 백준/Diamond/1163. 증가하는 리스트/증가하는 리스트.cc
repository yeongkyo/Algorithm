#include <iostream>
#include <string>
#include <vector>

using namespace std;

string smallest_match(string p) {
    if (p.empty()) return "";
    if (p[0] == '0') return "";
    string res = p;
    if (res[0] == '?') res[0] = '1';
    for (int i = 1; i < (int)res.length(); ++i) {
        if (res[i] == '?') res[i] = '0';
    }
    return res;
}

string get_smallest_greater(string p, string L) {
    if (p.length() < L.length()) return "";
    if (p.length() > L.length()) return smallest_match(p);
    
    string best = "";
    for (int idx = 0; idx < (int)p.length(); ++idx) {
        bool can_eq = true;
        string cand = p;
        for (int i = 0; i < idx; ++i) {
            if (p[i] == '?') cand[i] = L[i];
            else if (p[i] != L[i]) { can_eq = false; break; }
        }
        if (!can_eq) continue;
        
        if (p[idx] == '?') {
            if (L[idx] == '9') continue;
            cand[idx] = L[idx] + 1;
        } else {
            if (p[idx] <= L[idx]) continue;
        }
        
        for (int i = idx + 1; i < (int)p.length(); ++i) {
            if (cand[i] == '?') cand[i] = '0';
        }
        
        if (best == "" || cand < best) {
            best = cand;
        }
    }
    return best;
}

bool is_smaller_num_string(const string& a, const string& b) {
    if (a.length() != b.length()) return a.length() < b.length();
    return a < b;
}

bool can_form_valid(string S) {
    int n = S.length();
    vector<string> dp(n, "");
    
    for (int i = 0; i < n; ++i) {
        if (S[i] == ',') continue;
        
        for (int j = 0; j <= i; ++j) {
            bool has_comma = false;
            for(int k = j; k <= i; ++k) {
                if(S[k] == ',') {
                    has_comma = true;
                    break;
                }
            }
            if(has_comma) continue;
            
            if (j > 0 && S[j-1] != ',' && S[j-1] != '?') continue;
            
            string L = "0";
            if (j > 0) {
                if (j - 2 >= 0) {
                    L = dp[j-2];
                } else {
                    continue; 
                }
            }
            if (L == "") continue;
            
            string cand = get_smallest_greater(S.substr(j, i - j + 1), L);
            if (cand != "") {
                if (dp[i] == "" || is_smaller_num_string(cand, dp[i])) {
                    dp[i] = cand;
                }
            }
        }
    }
    return dp[n-1] != "";
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    string S_orig;
    if (!(cin >> S_orig)) return 0;
    
    string S = S_orig;
    int n = S.length();
    
    vector<char> chars = {',', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};
    
    for (int i = 0; i < n; ++i) {
        if (S_orig[i] != '?') continue;
        
        bool found = false;
        for (char c : chars) {
            if (c == ',' && i == 0) continue;
            if (c == '0' && (i == 0 || S[i-1] == ',')) continue;
            
            S[i] = c;
            if (can_form_valid(S)) {
                found = true;
                break;
            }
        }
        if (!found) {
            cout << -1 << "\n";
            return 0;
        }
    }
    
    if (can_form_valid(S)) {
        cout << S << "\n";
    } else {
        cout << -1 << "\n";
    }
    
    return 0;
}