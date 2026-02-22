#include <iostream>
#include <vector>

using namespace std;

int n;
int target[8];
int current_pts[8];
int remaining_matches[8];
vector<pair<int, int>> matches;

long long solve(int match_idx) {
    // 모든 경기를 다 치렀을 때
    if (match_idx == matches.size()) {
        for (int i = 0; i < n; ++i) {
            if (current_pts[i] != target[i]) return 0;
        }
        return 1;
    }

    int u = matches[match_idx].first;
    int v = matches[match_idx].second;

    // 현재 경기를 치르므로 남은 경기 수 1 감소
    remaining_matches[u]--;
    remaining_matches[v]--;

    long long ways = 0;

    // 1. u가 승리하는 경우 (u 3점, v 0점)
    if (current_pts[u] + 3 <= target[u] && 
        current_pts[u] + 3 + remaining_matches[u] * 3 >= target[u] &&
        current_pts[v] + remaining_matches[v] * 3 >= target[v]) {
        
        current_pts[u] += 3;
        ways += solve(match_idx + 1);
        current_pts[u] -= 3;
    }

    // 2. v가 승리하는 경우 (u 0점, v 3점)
    if (current_pts[v] + 3 <= target[v] && 
        current_pts[v] + 3 + remaining_matches[v] * 3 >= target[v] &&
        current_pts[u] + remaining_matches[u] * 3 >= target[u]) {
        
        current_pts[v] += 3;
        ways += solve(match_idx + 1);
        current_pts[v] -= 3;
    }

    // 3. 무승부인 경우 (u 1점, v 1점)
    if (current_pts[u] + 1 <= target[u] && current_pts[v] + 1 <= target[v] &&
        current_pts[u] + 1 + remaining_matches[u] * 3 >= target[u] &&
        current_pts[v] + 1 + remaining_matches[v] * 3 >= target[v]) {
        
        current_pts[u] += 1;
        current_pts[v] += 1;
        ways += solve(match_idx + 1);
        current_pts[u] -= 1;
        current_pts[v] -= 1;
    }

    // 백트래킹이 끝났으므로 남은 경기 수 원상복구
    remaining_matches[u]++;
    remaining_matches[v]++;

    return ways;
}

int main() {
    // 입출력 속도 향상
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    if (!(cin >> n)) return 0;
    
    for (int i = 0; i < n; ++i) {
        cin >> target[i];
        remaining_matches[i] = n - 1; // 초기에 각 팀은 (n-1) 경기를 치름
    }

    // 모든 팀들 간의 경기 대진표 생성
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            matches.push_back({i, j});
        }
    }

    cout << solve(0) << "\n";
    return 0;
}