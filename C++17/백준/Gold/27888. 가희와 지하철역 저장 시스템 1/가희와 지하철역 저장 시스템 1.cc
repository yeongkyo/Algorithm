#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <unordered_map>
#include <map>

using namespace std;

map<string, int> feature_to_id;
int feature_cnt = 0;
unordered_map<string, int> station_masks;
int mask_counts[512]; // 2^9 = 512

// 특징 문자열을 파싱하여 비트마스크로 변환
// is_query가 true일 때 모르는 특징이 나오면 -1 반환
int get_mask(string features, bool is_query) {
    int mask = 0;
    stringstream ss(features);
    string item;
    while (getline(ss, item, ',')) {
        if (feature_to_id.find(item) == feature_to_id.end()) {
            if (is_query) return -1; // 쿼리인데 모르는 특징이면 불가능한 조합
            feature_to_id[item] = feature_cnt++;
        }
        mask |= (1 << feature_to_id[item]);
    }
    return mask;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    if (!(cin >> n)) return 0;
    
    for (int i = 0; i < n; i++) {
        string name;
        cin >> name;
        station_masks[name] = 0;
        mask_counts[0]++;
    }

    int r;
    cin >> r;
    while (r--) {
        char cmd;
        cin >> cmd;
        if (cmd == 'U') {
            string sname, fstr;
            cin >> sname >> fstr;

            int old_mask = station_masks[sname];
            int new_mask = get_mask(fstr, false); // 새 특징 등록 가능

            mask_counts[old_mask]--;
            station_masks[sname] = new_mask;
            mask_counts[new_mask]++;
        } 
        else if (cmd == 'G') {
            string fstr;
            cin >> fstr;
            int target_mask = get_mask(fstr, true); // 모르는 특징 체크
            
            if (target_mask == -1) {
                cout << 0 << "\n";
                continue;
            }

            int total = 0;
            for (int i = 0; i < 512; i++) {
                if ((i & target_mask) == target_mask) {
                    total += mask_counts[i];
                }
            }
            cout << total << "\n";
        }
    }

    return 0;
}