#include <iostream>
#include <string>
#include <vector>

using namespace std;

int main() {
    string S;
    cin >> S;

    // 알파벳 26개(a-z)를 카운트할 배열을 0으로 초기화
    int counts[26] = {0, };

    // 단어 S를 한 글자씩 확인하며 해당 알파벳 위치의 카운트 증가
    for (char c : S) {
        counts[c - 'a']++;
    }

    // 결과 출력 (0번 인덱스부터 25번 인덱스까지 공백으로 구분)
    for (int i = 0; i < 26; i++) {
        cout << counts[i] << " ";
    }
    cout << endl;

    return 0;
}