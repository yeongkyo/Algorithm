#include <vector>
#include <algorithm>

using namespace std;

int n;
int min_total_cost;
vector<vector<int>> g_cost;
vector<vector<int>> g_hint;

void dfs(int stage_idx, int current_cost, vector<int>& coupon_counts) {
    if (current_cost >= min_total_cost) return;

    if (stage_idx == n - 1) {
        int used = min(coupon_counts[stage_idx], n - 1);
        int total = current_cost + g_cost[stage_idx][used];
        if (total < min_total_cost) {
            min_total_cost = total;
        }
        return;
    }

    int used = min(coupon_counts[stage_idx], n - 1);
    int stage_clear_cost = g_cost[stage_idx][used];

    dfs(stage_idx + 1, current_cost + stage_clear_cost, coupon_counts);

    int bundle_price = g_hint[stage_idx][0];

    for (size_t j = 1; j < g_hint[stage_idx].size(); ++j) {
        int target_stage = g_hint[stage_idx][j] - 1; 
        coupon_counts[target_stage]++;
    }

    dfs(stage_idx + 1, current_cost + stage_clear_cost + bundle_price, coupon_counts);

    for (size_t j = 1; j < g_hint[stage_idx].size(); ++j) {
        int target_stage = g_hint[stage_idx][j] - 1;
        coupon_counts[target_stage]--;
    }
}

int solution(vector<vector<int>> cost, vector<vector<int>> hint) {
    n = cost.size();
    g_cost = cost;
    g_hint = hint;
    min_total_cost = 2e9;

    vector<int> coupon_counts(n, 0);
    dfs(0, 0, coupon_counts);

    return min_total_cost;
}