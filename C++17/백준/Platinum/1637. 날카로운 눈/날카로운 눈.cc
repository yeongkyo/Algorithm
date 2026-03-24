#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

struct RangeInfo {
    long long a, c, b;
};

long long count_upto(const vector<RangeInfo>& ranges, long long x) {
    long long total = 0;
    for (const auto& r : ranges) {
        if (x < r.a) continue;
        long long end_point = min(r.c, x);
        total += (end_point - r.a) / r.b + 1;
    }
    return total;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    if (!(cin >> n)) return 0;

    vector<RangeInfo> ranges(n);
    for (int i = 0; i < n; i++) {
        cin >> ranges[i].a >> ranges[i].c >> ranges[i].b;
    }

    long long low = 1, high = 2147483647;
    
    if (count_upto(ranges, high) % 2 == 0) {
        cout << "NOTHING" << endl;
        return 0;
    }

    long long result_val = -1;
    while (low <= high) {
        long long mid = low + (high - low) / 2;
        if (count_upto(ranges, mid) % 2 != 0) {
            result_val = mid;
            high = mid - 1;
        } else {
            low = mid + 1;
        }
    }

    if (result_val == -1) {
        cout << "NOTHING" << endl;
    } else {
        long long final_count = count_upto(ranges, result_val) - count_upto(ranges, result_val - 1);
        cout << result_val << " " << final_count << endl;
    }

    return 0;
}