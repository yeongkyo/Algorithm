#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

int main(){
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    
    int N;
    cin >> N;
    
    vector<int> A(N), B(N);
    for (int i = 0; i < N; i++){
        cin >> A[i];
    }
    for (int i = 0; i < N; i++){
        cin >> B[i];
    }
    
    // B 수열에서 각 숫자가 등장하는 인덱스를 저장 (0-indexed)
    vector<int> pos(N + 1);
    for (int i = 0; i < N; i++){
        pos[B[i]] = i;
    }
    
    // A 수열을 B 수열에서의 인덱스로 변환
    vector<int> X(N);
    for (int i = 0; i < N; i++){
        X[i] = pos[A[i]];
    }
    
    // 변환된 수열 X의 최장 증가 부분 수열(LIS) 길이 계산
    vector<int> lis;
    for (int x : X) {
        auto it = lower_bound(lis.begin(), lis.end(), x);
        if(it == lis.end()){
            lis.push_back(x);
        } else {
            *it = x;
        }
    }
    
    cout << lis.size() << "\n";
    return 0;
}