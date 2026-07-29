#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>

// 주어진 상태에서 추가로 얻을 수 있는 최대 리프 노드 수를 반환하는 재귀 함수
static long long solve_dfs(long long cap, long long rem_dist, long long rem_split) {
    if (rem_dist <= 0 || rem_split < 2) {
        return 0;
    }
    
    long long max_additional_leaves = 0;
    
    // 선택지 1: 자식 수를 3개(s = 3)로 설정
    if (rem_split >= 3) {
        long long d = (cap < rem_dist) ? cap : rem_dist;
        long long leaves = d * 2;
        long long next_cap = d * 3;
        long long next_dist = rem_dist - d;
        long long next_split = rem_split / 3;
        
        long long future = solve_dfs(next_cap, next_dist, next_split);
        if (leaves + future > max_additional_leaves) {
            max_additional_leaves = leaves + future;
        }
    }
    
    // 선택지 2: 자식 수를 2개(s = 2)로 설정
    if (rem_split >= 2) {
        // 2a. 현재 가능한 최대로 배치
        long long d = (cap < rem_dist) ? cap : rem_dist;
        long long leaves = d * 1;
        long long next_cap = d * 2;
        long long next_dist = rem_dist - d;
        long long next_split = rem_split / 2;
        
        long long future = solve_dfs(next_cap, next_dist, next_split);
        if (leaves + future > max_additional_leaves) {
            max_additional_leaves = leaves + future;
        }
        
        // 2b. 다음 깊이에서 s = 3을 최대한 활용하기 위해 최소한만 배치
        if (cap >= rem_dist && rem_split >= 6) {
            long long d2 = (rem_dist + 2) / 3; // ceil(rem_dist / 3)
            if (d2 < d) {
                long long leaves2 = d2 * 1;
                long long next_cap2 = d2 * 2;
                long long next_dist2 = rem_dist - d2;
                long long next_split2 = rem_split / 2;
                
                long long future2 = solve_dfs(next_cap2, next_dist2, next_split2);
                if (leaves2 + future2 > max_additional_leaves) {
                    max_additional_leaves = leaves2 + future2;
                }
            }
        }
    }
    
    return max_additional_leaves;
}

int solution(int dist_limit, int split_limit) {
    if (dist_limit == 0) {
        return 1;
    }
    
    long long answer = 1 + solve_dfs(1, dist_limit, split_limit);
    return (int)answer;
}