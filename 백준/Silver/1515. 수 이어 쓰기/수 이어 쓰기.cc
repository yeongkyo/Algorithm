#include <iostream>
#include <string>
using namespace std;

// 주어진 N에 대해, "1 2 ... N"을 이어붙인 문자열 내에서 P가 부분 수열인지 확인하는 함수
bool canMatch(int N, const string &P) {
    int pIdx = 0, pSize = P.size();
    for (int i = 1; i <= N && pIdx < pSize; i++) {
        string num = to_string(i);
        for (char c : num) {
            if (c == P[pIdx]) {
                pIdx++;
                if (pIdx == pSize)
                    break;
            }
        }
    }
    return (pIdx == pSize);
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    
    string P;
    cin >> P;
    
    // 이분 탐색을 위한 범위 설정
    int low = 1, high = 2000000; // high는 충분히 큰 값으로 설정
    int ans = high;
    
    while(low <= high){
        int mid = low + (high - low) / 2;
        if(canMatch(mid, P)){
            ans = mid;
            high = mid - 1;
        } else {
            low = mid + 1;
        }
    }
    
    cout << ans;
    return 0;
}