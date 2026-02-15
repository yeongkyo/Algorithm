#include <bits/stdc++.h>
using namespace std;

struct Cow {
    int x, y;
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int N;
    cin >> N;
    vector<Cow> cows(N);
    vector<int> ys;
    ys.reserve(N);

    for (int i = 0; i < N; i++) {
        cin >> cows[i].x >> cows[i].y;
        ys.push_back(cows[i].y);
    }

    sort(cows.begin(), cows.end(), [](const Cow& a, const Cow& b){
        return a.x < b.x;
    });

    // y 후보는 y_i + 1만 보면 됨
    sort(ys.begin(), ys.end());
    ys.erase(unique(ys.begin(), ys.end()), ys.end());

    int best = N;

    for (int y0 : ys) {
        int b = y0 + 1; // fence y=b

        int leftBelow = 0, leftAbove = 0;
        int rightBelow = 0, rightAbove = 0;

        // 초기: 전부 right에 있음
        for (auto &c : cows) {
            if (c.y < b) rightBelow++;
            else rightAbove++;
        }

        // 세로선 스윕: 같은 x끼리 묶어서 right->left 이동
        int i = 0;
        while (i < N) {
            int xval = cows[i].x;

            // x = xval + 1 위치로 옮긴다고 생각하면,
            // x==xval 인 소들이 right->left로 넘어감
            int j = i;
            while (j < N && cows[j].x == xval) {
                if (cows[j].y < b) {
                    rightBelow--;
                    leftBelow++;
                } else {
                    rightAbove--;
                    leftAbove++;
                }
                j++;
            }

            int M = max({leftBelow, leftAbove, rightBelow, rightAbove});
            best = min(best, M);

            i = j;
        }
    }

    cout << best << "\n";
    return 0;
}
