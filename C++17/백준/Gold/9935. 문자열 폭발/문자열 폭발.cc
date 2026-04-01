#include <iostream>
#include <string>
#include <vector>

using namespace std;

int main() {
    // 입출력 속도 향상
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    string s, bomb;
    cin >> s >> bomb;

    string result = ""; // 스택 역할을 할 문자열
    int s_len = s.length();
    int b_len = bomb.length();

    for (int i = 0; i < s_len; i++) {
        result += s[i]; // 문자를 하나씩 추가

        // 현재 결과 문자열의 길이가 폭발 문자열보다 길거나 같을 때만 체크
        if (result.length() >= b_len) {
            bool match = true;
            
            // 끝에서부터 폭발 문자열과 일치하는지 확인
            for (int j = 0; j < b_len; j++) {
                if (result[result.length() - b_len + j] != bomb[j]) {
                    match = false;
                    break;
                }
            }

            // 일치하면 폭발 문자열 길이만큼 뒤에서 제거
            if (match) {
                for (int j = 0; j < b_len; j++) {
                    result.pop_back();
                }
            }
        }
    }

    // 결과 출력
    if (result.empty()) {
        cout << "FRULA" << "\n";
    } else {
        cout << result << "\n";
    }

    return 0;
}