#include <iostream>
#include <algorithm> // min 함수 사용을 위해 포함

using namespace std;

int main() {
    int burger1, burger2, burger3;
    int drink1, drink2;

    // 햄버거 3종 가격 입력
    cin >> burger1 >> burger2 >> burger3;
    // 음료 2종 가격 입력
    cin >> drink1 >> drink2;

    // 햄버거 중 최솟값 구하기
    int min_burger = min({burger1, burger2, burger3});
    
    // 음료 중 최솟값 구하기
    int min_drink = min(drink1, drink2);

    // 가장 싼 세트 메뉴 가격 계산 및 출력
    cout << min_burger + min_drink - 50 << endl;

    return 0;
}