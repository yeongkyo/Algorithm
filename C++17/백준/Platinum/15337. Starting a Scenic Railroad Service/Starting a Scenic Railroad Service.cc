#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

int board_cnt[100005];
int alight_cnt[100005];
int L[100005];
int R[100005];
int sweep[100005];

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int N;
    if (!(cin >> N)) return 0;

    vector<pair<int, int>> passengers(N);
    for (int i = 0; i < N; ++i) {
        cin >> passengers[i].first >> passengers[i].second;
        board_cnt[passengers[i].first]++;
        alight_cnt[passengers[i].second]++;
        sweep[passengers[i].first]++;
        sweep[passengers[i].second]--;
    }

    int current_alight = 0;
    for (int i = 1; i <= 100000; ++i) {
        current_alight += alight_cnt[i];
        L[i] = current_alight;
    }

    int current_board = 0;
    for (int i = 100000; i >= 1; --i) {
        current_board += board_cnt[i];
        R[i] = current_board;
    }

    int s1 = 0;
    for (int i = 0; i < N; ++i) {
        s1 = max(s1, N - L[passengers[i].first] - R[passengers[i].second]);
    }

    int s2 = 0;
    int current_passengers = 0;
    for (int i = 1; i <= 100000; ++i) {
        current_passengers += sweep[i];
        s2 = max(s2, current_passengers);
    }

    cout << s1 << " " << s2 << endl;

    return 0;
}