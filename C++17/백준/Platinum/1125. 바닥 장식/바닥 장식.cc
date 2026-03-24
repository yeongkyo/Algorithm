#include <iostream>
#include <algorithm>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    long long x1, y1, x2, y2;
    if (!(cin >> x1 >> y1 >> x2 >> y2)) return 0;

    long long C_even[6] = {0}, C_odd[6] = {0};
    for (long long i = x1 / 5; i <= (x2 - 1) / 5; ++i) {
        long long w = min(x2, i * 5 + 5) - max(x1, i * 5);
        if (w > 0) {
            if (i % 2 == 0) C_even[w]++;
            else C_odd[w]++;
        }
    }

    long long R_even[6] = {0}, R_odd[6] = {0};
    for (long long j = y1 / 5; j <= (y2 - 1) / 5; ++j) {
        long long h = min(y2, j * 5 + 5) - max(y1, j * 5);
        if (h > 0) {
            if (j % 2 == 0) R_even[h]++;
            else R_odd[h]++;
        }
    }

    long long counts[6] = {0};
    for (int w = 1; w <= 5; ++w) {
        for (int h = 1; h <= 5; ++h) {
            counts[w] += C_even[w] * R_even[h] * h;
            counts[w] += C_odd[w] * R_odd[h] * h;
            counts[h] += C_even[w] * R_odd[h] * w;
            counts[h] += C_odd[w] * R_even[h] * w;
        }
    }

    long long bins = 0;
    
    bins += counts[5];

    bins += counts[4];
    counts[1] = max(0LL, counts[1] - counts[4]);

    long long p32 = min(counts[3], counts[2]);
    bins += p32;
    counts[3] -= p32;
    counts[2] -= p32;

    if (counts[3] > 0) {
        bins += counts[3];
        counts[1] = max(0LL, counts[1] - counts[3] * 2);
    }

    if (counts[2] > 0) {
        long long b2 = counts[2] / 2;
        bins += b2;
        // 길이가 2인 조각 2개를 하나의 5짜리 판에 넣으면 1칸이 남음
        counts[1] = max(0LL, counts[1] - b2);

        if (counts[2] % 2 != 0) {
            bins += 1;
            counts[1] = max(0LL, counts[1] - 3);
        }
    }

    if (counts[1] > 0) {
        bins += (counts[1] + 4) / 5;
    }

    cout << bins << "\n";

    return 0;
}