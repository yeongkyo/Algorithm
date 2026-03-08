#include <iostream>
#include <vector>
#include <numeric>

int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(NULL);
    std::cout.tie(NULL);

    int n, x;
    std::cin >> n >> x;

    std::vector<int> visitors(n);
    for (int i = 0; i < n; ++i) {
        std::cin >> visitors[i];
    }

    long long current_sum = 0;
    for (int i = 0; i < x; ++i) {
        current_sum += visitors[i];
    }

    long long max_sum = current_sum;
    int count = 1;

    if (n == x) { // N과 X가 같을 경우, 슬라이딩 윈도우 로직을 탈 필요 없음
        if (max_sum == 0) {
            std::cout << "SAD\n";
        } else {
            std::cout << max_sum << "\n";
            std::cout << count << "\n";
        }
        return 0;
    }
    
    for (int i = x; i < n; ++i) {
        current_sum += visitors[i];
        current_sum -= visitors[i - x];

        if (current_sum > max_sum) {
            max_sum = current_sum;
            count = 1;
        } else if (current_sum == max_sum) {
            count++;
        }
    }

    if (max_sum == 0) {
        std::cout << "SAD\n";
    } else {
        std::cout << max_sum << "\n";
        std::cout << count << "\n";
    }

    return 0;
}