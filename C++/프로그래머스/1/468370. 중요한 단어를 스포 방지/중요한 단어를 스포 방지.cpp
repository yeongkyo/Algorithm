#include <string>
#include <vector>
#include <unordered_set>
#include <algorithm>

using namespace std;

struct WordInfo {
    string text;
    int start;
    int end;
    int max_r; 
};

int solution(string message, vector<vector<int>> spoiler_ranges) {
    int n = message.length();
    int m = spoiler_ranges.size();
    
    vector<WordInfo> words;
    int word_start = -1;
    
    for (int i = 0; i < n; ++i) {
        if (message[i] != ' ') {
            if (word_start == -1) word_start = i;
        } else {
            if (word_start != -1) {
                words.push_back({message.substr(word_start, i - word_start), word_start, i - 1, -1});
                word_start = -1;
            }
        }
    }
    if (word_start != -1) {
        words.push_back({message.substr(word_start, n - word_start), word_start, n - 1, -1});
    }

    unordered_set<string> non_spoiler_words;
    vector<vector<int>> range_to_words(m); 
    
    int r_idx = 0;
    for (int w = 0; w < (int)words.size(); ++w) {
        int w_start = words[w].start;
        int w_end = words[w].end;
        
        while (r_idx < m && spoiler_ranges[r_idx][1] < w_start) {
            r_idx++;
        }
        
        int max_r = -1;
        for (int r = r_idx; r < m; ++r) {
            if (spoiler_ranges[r][0] > w_end) break; 
            
            if (max(w_start, spoiler_ranges[r][0]) <= min(w_end, spoiler_ranges[r][1])) {
                max_r = max(max_r, r);
            }
        }

        words[w].max_r = max_r;
        
        if (max_r == -1) {
            non_spoiler_words.insert(words[w].text);
        } else {
            range_to_words[max_r].push_back(w);
        }
    }

    unordered_set<string> revealed_spoiler_words;
    int answer = 0;

    for (int r = 0; r < m; ++r) {
        for (int w_idx : range_to_words[r]) {
            const string& word = words[w_idx].text;
            
            if (non_spoiler_words.count(word)) continue;
            
            if (revealed_spoiler_words.count(word)) continue;

            answer++;
            revealed_spoiler_words.insert(word);
        }
    }

    return answer;
}