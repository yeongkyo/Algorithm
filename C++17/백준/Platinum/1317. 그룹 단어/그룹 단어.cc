#include <iostream>
#include <vector>
#include <string>

using namespace std;

bool isValid(const string& s) {
    vector<bool> seen(26, false);
    for (int i = 0; i < s.length(); ++i) {
        if (i > 0 && s[i] != s[i - 1]) {
            if (seen[s[i] - 'a']) return false;
        }
        seen[s[i] - 'a'] = true;
    }
    return true;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    if (!(cin >> n)) return 0;

    vector<string> pieces(n);
    for (int i = 0; i < n; ++i) {
        cin >> pieces[i];
        if (!isValid(pieces[i])) {
            cout << "gg\n";
            return 0;
        }
    }

    bool merged_any = true;
    while (merged_any) {
        merged_any = false;
        for (char c = 'a'; c <= 'z'; ++c) {
            int L_idx = -1;
            int R_idx = -1;
            vector<int> pure_idx;
            int internal_count = 0;
            int total_c_pieces = 0;

            for (int i = 0; i < pieces.size(); ++i) {
                if (pieces[i].find(c) == string::npos) continue;

                total_c_pieces++;
                
                bool all_c = true;
                for (char ch : pieces[i]) {
                    if (ch != c) { all_c = false; break; }
                }

                if (all_c) {
                    pure_idx.push_back(i);
                } else if (pieces[i].front() == c) {
                    if (L_idx != -1) { cout << "gg\n"; return 0; }
                    L_idx = i;
                } else if (pieces[i].back() == c) {
                    if (R_idx != -1) { cout << "gg\n"; return 0; }
                    R_idx = i;
                } else {
                    internal_count++;
                }
            }

            if (total_c_pieces <= 1) continue;

            if (internal_count > 0) {
                cout << "gg\n";
                return 0;
            }

            string merged = "";
            if (R_idx != -1) merged += pieces[R_idx];
            for (int idx : pure_idx) merged += pieces[idx];
            if (L_idx != -1) merged += pieces[L_idx];

            if (!isValid(merged)) {
                cout << "gg\n";
                return 0;
            }

            vector<string> next_pieces;
            for (int i = 0; i < pieces.size(); ++i) {
                if (i == R_idx || i == L_idx) continue;
                bool is_pure = false;
                for (int p_idx : pure_idx) {
                    if (i == p_idx) { is_pure = true; break; }
                }
                if (is_pure) continue;
                next_pieces.push_back(pieces[i]);
            }
            
            next_pieces.push_back(merged);
            pieces = next_pieces;
            merged_any = true;
            break;
        }
    }

    if (pieces.size() > 1) {
        cout << "-_-\n";
    } else {
        cout << pieces[0] << "\n";
    }

    return 0;
}