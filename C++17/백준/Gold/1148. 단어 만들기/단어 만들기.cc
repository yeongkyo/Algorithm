#include <iostream>
#include <vector>
#include <string>

using namespace std;

struct Word {
    int mask;
    int8_t cnt[26];
};

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    vector<Word> dict;
    string s;
    
    while (cin >> s && s != "-") {
        Word w = {0, {0}};
        for (char c : s) {
            w.cnt[c - 'A']++;
            w.mask |= (1 << (c - 'A'));
        }
        dict.push_back(w);
    }

    while (cin >> s && s != "#") {
        int8_t board_cnt[26] = {0};
        bool present[26] = {false};
        int board_mask = 0;

        for (char c : s) {
            board_cnt[c - 'A']++;
            present[c - 'A'] = true;
            board_mask |= (1 << (c - 'A'));
        }

        int valid_counts[26] = {0};

        for (const auto& w : dict) {
            if ((w.mask & ~board_mask) != 0) continue;

            bool ok = true;
            for (int i = 0; i < 26; ++i) {
                if (w.cnt[i] > board_cnt[i]) {
                    ok = false;
                    break;
                }
            }

            if (ok) {
                for (int i = 0; i < 26; ++i) {
                    if (w.cnt[i] > 0) {
                        valid_counts[i]++;
                    }
                }
            }
        }

        int min_cnt = 2e9, max_cnt = -1;
        for (int i = 0; i < 26; ++i) {
            if (present[i]) {
                if (valid_counts[i] < min_cnt) min_cnt = valid_counts[i];
                if (valid_counts[i] > max_cnt) max_cnt = valid_counts[i];
            }
        }

        string min_chars = "", max_chars = "";
        for (int i = 0; i < 26; ++i) {
            if (present[i]) {
                if (valid_counts[i] == min_cnt) min_chars += (char)(i + 'A');
                if (valid_counts[i] == max_cnt) max_chars += (char)(i + 'A');
            }
        }

        cout << min_chars << " " << min_cnt << " " << max_chars << " " << max_cnt << "\n";
    }

    return 0;
}