#include <string>
#include <vector>
#include <algorithm>

using namespace std;

extern string submit(int);

struct Candidate {
    int num;
    int d[4];
};

pair<int, int> calc_feedback(const int g[4], const int c[4]) {
    int S = 0, B = 0;
    for (int i = 0; i < 4; i++) {
        if (g[i] == c[i]) {
            S++;
        } else {
            for (int j = 0; j < 4; j++) {
                if (g[i] == c[j]) {
                    B++;
                    break;
                }
            }
        }
    }
    return {S, B};
}

int solution(int n) {
    vector<Candidate> all_cands;
    for (int a = 1; a <= 9; a++) {
        for (int b = 1; b <= 9; b++) {
            if (b == a) continue;
            for (int c = 1; c <= 9; c++) {
                if (c == a || c == b) continue;
                for (int d = 1; d <= 9; d++) {
                    if (d == a || d == b || d == c) continue;
                    all_cands.push_back({a * 1000 + b * 100 + c * 10 + d, {a, b, c, d}});
                }
            }
        }
    }

    vector<Candidate> C = all_cands;

    while (!C.empty()) {
        int guess_num = 0;
        int guess_d[4];

        if (C.size() == all_cands.size()) {
            guess_num = 1234;
            guess_d[0] = 1; guess_d[1] = 2; guess_d[2] = 3; guess_d[3] = 4;
        } else if (C.size() == 1) {
            guess_num = C[0].num;
            for (int i = 0; i < 4; i++) guess_d[i] = C[0].d[i];
        } else {
            int best_max_group = 1e9;
            int best_sum_sq = 1e9;
            bool best_in_C = false;

            for (const auto& g : all_cands) {
                int counts[5][5] = {0};
                for (const auto& c : C) {
                    auto [S, B] = calc_feedback(g.d, c.d);
                    counts[S][B]++;
                }

                int max_g = 0;
                int sum_sq = 0;
                for (int s = 0; s <= 4; s++) {
                    for (int b = 0; b <= 4; b++) {
                        max_g = max(max_g, counts[s][b]);
                        sum_sq += counts[s][b] * counts[s][b];
                    }
                }

                bool in_C = false;
                for (const auto& c : C) {
                    if (c.num == g.num) {
                        in_C = true;
                        break;
                    }
                }

                bool is_better = false;
                if (max_g < best_max_group) {
                    is_better = true;
                } else if (max_g == best_max_group) {
                    if (in_C && !best_in_C) {
                        is_better = true;
                    } else if (in_C == best_in_C) {
                        if (sum_sq < best_sum_sq) {
                            is_better = true;
                        }
                    }
                }

                if (is_better) {
                    best_max_group = max_g;
                    best_sum_sq = sum_sq;
                    best_in_C = in_C;
                    guess_num = g.num;
                    for (int i = 0; i < 4; i++) guess_d[i] = g.d[i];
                }
            }
        }

        string res = submit(guess_num);
        if (res == "4S 0B") {
            return guess_num;
        }

        int res_S = res[0] - '0';
        int res_B = res[3] - '0';

        vector<Candidate> next_C;
        for (const auto& c : C) {
            auto [S, B] = calc_feedback(guess_d, c.d);
            if (S == res_S && B == res_B) {
                next_C.push_back(c);
            }
        }
        C = next_C;
    }

    return 0;
}