#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int N, M;
    if (!(cin >> N >> M)) return 0;

    vector<vector<int>> cards(N, vector<int>(M));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < M; ++j) {
            cin >> cards[i][j];
        }
    }

    int min_total_moves = 1e9;

    for (int joker = 0; joker < N; ++joker) {
        int non_empty_count = 0;
        vector<bool> color_has_pure_box(M, false);

        for (int i = 0; i < N; ++i) {
            if (i == joker) continue; 

            int total_in_box = 0;
            int distinct_colors = 0;
            int last_color = -1;

            for (int j = 0; j < M; ++j) {
                if (cards[i][j] > 0) {
                    total_in_box += cards[i][j];
                    distinct_colors++;
                    last_color = j;
                }
            }

            if (total_in_box > 0) {
                non_empty_count++;
                if (distinct_colors == 1) {
                    color_has_pure_box[last_color] = true;
                }
            }
        }

        int pure_colors_count = 0;
        for (int j = 0; j < M; ++j) {
            if (color_has_pure_box[j]) pure_colors_count++;
        }

        int current_moves = non_empty_count - pure_colors_count;
        if (current_moves < min_total_moves) {
            min_total_moves = current_moves;
        }
    }

    if (min_total_moves == 1e9) cout << -1 << endl;
    else cout << min_total_moves << endl;

    return 0;
}