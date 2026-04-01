#include <iostream>
#include <stack>
#include <vector>

using namespace std;

// 탑의 정보를 저장할 구조체 (높이와 번호)
struct Tower {
    int height;
    int index;
};

int main() {
    // 입출력 속도 향상
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    int N;
    cin >> N;

    stack<Tower> s;

    for (int i = 1; i <= N; i++) {
        int h;
        cin >> h;

        // 1. 현재 탑보다 낮은 앞선 탑들은 스택에서 제거
        // 어차피 현재 탑에 가려져서 의미가 없어짐
        while (!s.empty() && s.top().height < h) {
            s.pop();
        }

        // 2. 수신 가능한 탑 확인
        if (s.empty()) {
            cout << 0 << " ";
        } else {
            // 스택의 top에 있는 탑이 나보다 높으면서 가장 가까운 탑
            cout << s.top().index << " ";
        }

        // 3. 현재 탑을 스택에 추가
        s.push({h, i});
    }

    cout << "\n";

    return 0;
}