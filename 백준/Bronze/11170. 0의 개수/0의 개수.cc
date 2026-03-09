#include <iostream>

using namespace std;

// 숫자에 포함된 0의 개수를 구하는 함수
long long countZeros(int n) {
    if (n == 0) return 1;
    
    long long count = 0;
    while (n > 0) {
        if (n % 10 == 0) {
            count++;
        }
        n /= 10;
    }
    return count;
}

int main() {
    // 입출력 속도 향상
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int T;
    cin >> T;

    while (T--) {
        int N, M;
        cin >> N >> M;

        long long totalZeros = 0;
        for (int i = N; i <= M; ++i) {
            totalZeros += countZeros(i);
        }
        cout << totalZeros << "\n";
    }

    return 0;
}