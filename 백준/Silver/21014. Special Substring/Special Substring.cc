#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    int n, k;
    if (!(cin >> n >> k)) return 0;

    string s;
    cin >> s;

    int max_freq_overall = 0;

    for (char target = 'A'; target <= 'Z'; ++target) {
        int current_window_count = 0;

        for (int i = 0; i < k; ++i) {
            if (s[i] == target) {
                current_window_count++;
            }
        }
        max_freq_overall = max(max_freq_overall, current_window_count);

        for (int i = k; i < n; ++i) {
            if (s[i] == target) current_window_count++;
            if (s[i - k] == target) current_window_count--;
            
            max_freq_overall = max(max_freq_overall, current_window_count);
        }
        
        if (max_freq_overall == k) break;
    }

    cout << k - max_freq_overall << "\n";

    return 0;
}