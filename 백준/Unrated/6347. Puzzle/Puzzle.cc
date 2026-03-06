#include <iostream>
#include <vector>
#include <string>

using namespace std;

int piece_count[81];
int grid_R[36];
int grid_B[36];
int N, M;

int char_to_int(char c) {
    if (c == 'F') return 0;
    if (c == 'O') return 1;
    if (c == 'I') return 2;
    return 0;
}

bool solve(int idx) {
    if (idx == N * M) return true; 

    int r = idx / M;
    int c = idx % M;

    int req_T = (r == 0) ? 0 : 3 - grid_B[idx - M];
    int req_L = (c == 0) ? 0 : 3 - grid_R[idx - 1];

    for (int R = 0; R <= 2; ++R) {
        if ((c == M - 1) && R != 0) continue; 
        if ((c < M - 1) && R == 0) continue; 
        for (int B = 0; B <= 2; ++B) {
            if ((r == N - 1) && B != 0) continue; 
            if ((r < N - 1) && B == 0) continue; 

            int p_idx = req_T * 27 + R * 9 + B * 3 + req_L;
            
            if (piece_count[p_idx] > 0) {
                piece_count[p_idx]--;
                grid_R[idx] = R;
                grid_B[idx] = B;
                
                if (solve(idx + 1)) return true;
                
                piece_count[p_idx]++;
            }
        }
    }
    return false;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    while (cin >> N >> M && (N != 0 || M != 0)) {
        for (int i = 0; i < 81; ++i) piece_count[i] = 0;
        
        int f_top = 0, f_bottom = 0, f_left = 0, f_right = 0;
        int o_right = 0, i_left = 0;
        int i_right = 0, o_left = 0;
        int o_bottom = 0, i_top = 0;
        int i_bottom = 0, o_top = 0;

        for (int i = 0; i < N * M; ++i) {
            string s;
            cin >> s;
            int T = char_to_int(s[0]);
            int R = char_to_int(s[1]);
            int B = char_to_int(s[2]);
            int L = char_to_int(s[3]);
            
            piece_count[T * 27 + R * 9 + B * 3 + L]++;
            
            if (T == 0) f_top++; else if (T == 1) o_top++; else i_top++;
            if (R == 0) f_right++; else if (R == 1) o_right++; else i_right++;
            if (B == 0) f_bottom++; else if (B == 1) o_bottom++; else i_bottom++;
            if (L == 0) f_left++; else if (L == 1) o_left++; else i_left++;
        }
        
        bool possible = true;
        if (f_top != M || f_bottom != M || f_left != N || f_right != N) possible = false;
        if (o_right != i_left || i_right != o_left) possible = false;
        if (o_bottom != i_top || i_bottom != o_top) possible = false;
        
        if (possible) {
            if (solve(0)) cout << "YES\n";
            else cout << "NO\n";
        } else {
            cout << "NO\n";
        }
    }
    return 0;
}