#include <iostream>

using namespace std;

int main() {
    int a, b, c;
    
    // 세 줄에 걸쳐 각을 입력받음
    cin >> a >> b >> c;

    // 1. 세 각의 합이 180이 아닌 경우
    if (a + b + c != 180) {
        cout << "Error" << endl;
    } 
    // 2. 세 각의 합이 180인 경우 중 세 각이 모두 60인 경우
    else if (a == 60 && b == 60 && c == 60) {
        cout << "Equilateral" << endl;
    }
    // 3. 세 각의 합이 180인 경우 중 두 각이 같은 경우
    // (a==b 또는 b==c 또는 a==c)
    else if (a == b || b == c || a == c) {
        cout << "Isosceles" << endl;
    }
    // 4. 세 각의 합이 180인 경우 중 같은 각이 없는 경우
    else {
        cout << "Scalene" << endl;
    }

    return 0;
}