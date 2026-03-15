#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>

using namespace std;

int R_pts[2500][2], G_pts[2500][2], B_pts[2500][2];
int num_R = 0, num_G = 0, num_B = 0;

int min_f[3][105][105];
int max_f[3][105][105];
vector<int> ext_pts[3][105][105];

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);
    
    int N, M;
    if (!(cin >> N >> M)) return 0;
    
    for (int r = 0; r < N; ++r) {
        string row;
        cin >> row;
        for (int c = 0; c < M; ++c) {
            if (row[c] == 'R') {
                R_pts[num_R][0] = r;
                R_pts[num_R][1] = c;
                num_R++;
            } else if (row[c] == 'G') {
                G_pts[num_G][0] = r;
                G_pts[num_G][1] = c;
                num_G++;
            } else if (row[c] == 'B') {
                B_pts[num_B][0] = r;
                B_pts[num_B][1] = c;
                num_B++;
            }
        }
    }
    
    if (num_R == 0 || num_G == 0 || num_B == 0) {
        cout << 0 << "\n";
        return 0;
    }
    
    for (int color = 0; color < 3; ++color) {
        int num_pts = (color == 0) ? num_R : ((color == 1) ? num_G : num_B);
        int (*pts)[2] = (color == 0) ? R_pts : ((color == 1) ? G_pts : B_pts);
        
        for (int nc = -50; nc <= 50; ++nc) {
            for (int nr = -50; nr <= 50; ++nr) {
                if (nc == 0 && nr == 0) continue;
                int mn = 1e9, mx = -1e9;
                for (int i = 0; i < num_pts; ++i) {
                    int val = nc * pts[i][1] + nr * pts[i][0];
                    if (val < mn) mn = val;
                    if (val > mx) mx = val;
                }
                min_f[color][nc+50][nr+50] = mn;
                max_f[color][nc+50][nr+50] = mx;
                
                for (int i = 0; i < num_pts; ++i) {
                    int val = nc * pts[i][1] + nr * pts[i][0];
                    if (val == mn || val == mx) {
                        ext_pts[color][nc+50][nr+50].push_back(i);
                    }
                }
            }
        }
    }
    
    long long non_enlargeable = 0;
    for (int i = 0; i < num_R; ++i) {
        int r1 = R_pts[i][0];
        int c1 = R_pts[i][1];
        for (int j = 0; j < num_G; ++j) {
            int r2 = G_pts[j][0];
            int c2 = G_pts[j][1];
            
            int nc1 = r1 - r2;
            int nr1 = c2 - c1;
            
            int C1 = nc1 * c1 + nr1 * r1;
            
            int mn1 = min_f[2][nc1+50][nr1+50];
            int mx1 = max_f[2][nc1+50][nr1+50];
            int max_val1 = max(abs(mn1 - C1), abs(mx1 - C1));
            
            for (int k : ext_pts[2][nc1+50][nr1+50]) {
                int r3 = B_pts[k][0];
                int c3 = B_pts[k][1];
                
                int val1 = nc1 * c3 + nr1 * r3;
                if (abs(val1 - C1) != max_val1) continue;
                
                int nc2 = r2 - r3;
                int nr2 = c3 - c2;
                int C2 = nc2 * c2 + nr2 * r2;
                int mn2 = min_f[0][nc2+50][nr2+50];
                int mx2 = max_f[0][nc2+50][nr2+50];
                int max_val2 = max(abs(mn2 - C2), abs(mx2 - C2));
                if (max_val1 != max_val2) continue;
                
                int nc3 = r3 - r1;
                int nr3 = c1 - c3;
                int C3 = nc3 * c3 + nr3 * r3;
                int mn3 = min_f[1][nc3+50][nr3+50];
                int mx3 = max_f[1][nc3+50][nr3+50];
                int max_val3 = max(abs(mn3 - C3), abs(mx3 - C3));
                if (max_val1 != max_val3) continue;
                
                non_enlargeable++;
            }
        }
    }
    
    long long total_beautiful = (long long)num_R * num_G * num_B;
    long long ans = total_beautiful - non_enlargeable;
    cout << ans << "\n";
    
    return 0;
}