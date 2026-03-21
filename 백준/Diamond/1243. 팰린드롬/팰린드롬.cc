#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <cstring>

using namespace std;

int N, L;
vector<string> words;
vector<string> states;
unordered_map<string, int> str_id;

int get_id(const string& s) {
    auto it = str_id.find(s);
    if (it == str_id.end()) {
        str_id[s] = states.size();
        states.push_back(s);
        return states.size() - 1;
    }
    return it->second;
}

struct Trans {
    bool valid;
    int nxt_id;
    int nxt_side;
};

vector<vector<vector<Trans>>> nxt;

bool is_palindrome(const string& s) {
    int i = 0, j = (int)s.length() - 1;
    while (i < j) {
        if (s[i] != s[j]) return false;
        i++; j--;
    }
    return true;
}

long long dp[35][1600][2];

long long solve(int len, int id, int side) {
    if (len == L) {
        return is_palindrome(states[id]) ? 1 : 0;
    }
    if (len > L) return 0;
    if (dp[len][id][side] != -1) return dp[len][id][side];
    
    long long ans = 0;
    for (int i = 0; i < N; ++i) {
        if (nxt[id][i][side].valid) {
            ans += solve(len + words[i].length(), nxt[id][i][side].nxt_id, nxt[id][i][side].nxt_side);
        }
    }
    return dp[len][id][side] = ans;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    if (!(cin >> N >> L)) return 0;
    
    words.resize(N);
    for (int i = 0; i < N; ++i) {
        cin >> words[i];
    }
    
    get_id("");
    for (const string& w : words) {
        for (int i = 0; i < w.length(); ++i) {
            get_id(w.substr(i));
        }
        string wr = w;
        reverse(wr.begin(), wr.end());
        for (int i = 0; i < wr.length(); ++i) {
            get_id(wr.substr(i));
        }
    }
    
    int num_states = states.size();
    nxt.assign(num_states, vector<vector<Trans>>(N, vector<Trans>(2)));
    
    for (int id = 0; id < num_states; ++id) {
        string req = states[id];
        for (int i = 0; i < N; ++i) {
            string w = words[i];
            
            // side == 1 (왼쪽에 단어 추가)
            int M = min(w.length(), req.length());
            if (w.substr(0, M) == req.substr(0, M)) {
                nxt[id][i][1].valid = true;
                if (w.length() < req.length()) {
                    nxt[id][i][1].nxt_id = get_id(req.substr(M));
                    nxt[id][i][1].nxt_side = 1;
                } else if (w.length() > req.length()) {
                    nxt[id][i][1].nxt_id = get_id(w.substr(M));
                    nxt[id][i][1].nxt_side = 0;
                } else {
                    nxt[id][i][1].nxt_id = get_id("");
                    nxt[id][i][1].nxt_side = 1;
                }
            } else {
                nxt[id][i][1].valid = false;
            }
            
            // side == 0 (오른쪽에 단어 추가)
            string wr = w;
            reverse(wr.begin(), wr.end());
            M = min(wr.length(), req.length());
            if (wr.substr(0, M) == req.substr(0, M)) {
                nxt[id][i][0].valid = true;
                if (wr.length() < req.length()) {
                    nxt[id][i][0].nxt_id = get_id(req.substr(M));
                    nxt[id][i][0].nxt_side = 0;
                } else if (wr.length() > req.length()) {
                    nxt[id][i][0].nxt_id = get_id(wr.substr(M));
                    nxt[id][i][0].nxt_side = 1;
                } else {
                    nxt[id][i][0].nxt_id = get_id("");
                    nxt[id][i][0].nxt_side = 1;
                }
            } else {
                nxt[id][i][0].valid = false;
            }
        }
    }
    
    memset(dp, -1, sizeof(dp));
    
    long long result = solve(0, get_id(""), 1);
    cout << result << "\n";
    
    return 0;
}