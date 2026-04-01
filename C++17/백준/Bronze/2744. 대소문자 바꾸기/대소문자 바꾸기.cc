#include <iostream>
#include <string>
#include <cctype> // isupper, islower, toupper, tolower 함수 사용을 위해 포함

using namespace std;

int main() {
    string s;
    cin >> s;

    for (int i = 0; i < s.length(); i++) {
        // 대문자라면 소문자로 변환
        if (isupper(s[i])) {
            s[i] = tolower(s[i]);
        }
        // 소문자라면 대문자로 변환
        else if (islower(s[i])) {
            s[i] = toupper(s[i]);
        }
    }

    cout << s << endl;

    return 0;
}