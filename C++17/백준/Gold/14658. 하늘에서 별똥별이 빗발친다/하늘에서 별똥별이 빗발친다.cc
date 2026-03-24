#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

struct Point {
    int x, y;
};

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int N, M, L, K;
    cin >> N >> M >> L >> K;

    vector<Point> stars(K);
    for (int i = 0; i < K; ++i) {
        cin >> stars[i].x >> stars[i].y;
    }

    int max_bounced = 0;

    for (int i = 0; i < K; ++i) {
        for (int j = 0; j < K; ++j) {
            int start_x = stars[i].x;
            int start_y = stars[j].y;
            int current_bounced = 0;

            for (int k = 0; k < K; ++k) {
                if (start_x <= stars[k].x && stars[k].x <= start_x + L &&
                    start_y <= stars[k].y && stars[k].y <= start_y + L) {
                    current_bounced++;
                }
            }

            if (current_bounced > max_bounced) {
                max_bounced = current_bounced;
            }
        }
    }

    cout << K - max_bounced << "\n";

    return 0;
}