#include <iostream>
#include <vector>
#include <string>

using namespace std;

int N;
vector<string> board;

int dy[] = {-1, -1, -1, 0, 0, 1, 1, 1};
int dx[] = {-1, 0, 1, -1, 1, -1, 0, 1};

int countFlips(int y, int x) {
    int total = 0;

    for (int i = 0; i < 8; i++) {
        int ny = y + dy[i];
        int nx = x + dx[i];
        int count = 0;

        while (ny >= 0 && ny < N && nx >= 0 && nx < N && board[ny][nx] == 'W') {
            ny += dy[i];
            nx += dx[i];
            count++;
        }

        if (ny >= 0 && ny < N && nx >= 0 && nx < N && board[ny][nx] == 'B') {
            total += count;
        }
    }
    return total;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    if (!(cin >> N)) return 0;

    board.resize(N);
    for (int i = 0; i < N; i++) {
        cin >> board[i];
    }

    int maxFlips = 0;
    int bestX = -1, bestY = -1;

    for (int y = 0; y < N; y++) {
        for (int x = 0; x < N; x++) {
            if (board[y][x] == '.') {
                int flipped = countFlips(y, x);
                
                if (flipped > maxFlips) {
                    maxFlips = flipped;
                    bestX = x;
                    bestY = y;
                }
            }
        }
    }

    if (maxFlips == 0) {
        cout << "PASS" << endl;
    } else {
        cout << bestX << " " << bestY << endl;
        cout << maxFlips << endl;
    }

    return 0;
}