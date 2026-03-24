#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

struct Node {
    long long weight;
    int y;
};

int n;
long long weights[262145];
Node nodes[262145];
int x_idx = 1;

void inorder(int node, int depth) {
    if (node > n) return;
    
    inorder(node * 2, depth + 1);
    nodes[x_idx] = {weights[node], depth};
    x_idx++;
    inorder(node * 2 + 1, depth + 1);
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    if (!(cin >> n)) return 0;

    long long max_single_weight = -1e18;
    for (int i = 1; i <= n; i++) {
        cin >> weights[i];
        max_single_weight = max(max_single_weight, weights[i]);
    }

    if (max_single_weight <= 0) {
        cout << max_single_weight << "\n";
        return 0;
    }

    inorder(1, 0);

    int max_depth = 0;
    for(int i = 1; i <= n; i++) {
        max_depth = max(max_depth, nodes[i].y);
    }

    long long max_sum = -1e18;

    for (int top = 0; top <= max_depth; top++) {
        for (int bottom = top; bottom <= max_depth; bottom++) {
            
            long long current_sum = 0;
            for (int i = 1; i <= n; i++) {
                if (nodes[i].y >= top && nodes[i].y <= bottom) {
                    current_sum += nodes[i].weight;
                    if (current_sum < 0) current_sum = 0;
                    max_sum = max(max_sum, current_sum);
                }
            }
            
        }
    }

    cout << max_sum << "\n";
    return 0;
}