#include <iostream>
#include <vector>
#include <numeric>

using namespace std;

// 부모 노드를 저장할 배열
int parent[1000001];

// Find 연산 (경로 압축 적용)
int find(int x) {
    if (parent[x] == x) return x;
    return parent[x] = find(parent[x]); // 부모를 루트로 바로 연결
}

// Union 연산
void unite(int a, int b) {
    a = find(a);
    b = find(b);
    if (a != b) {
        parent[a] = b; // a의 루트를 b의 루트 아래로 병합
    }
}

int main() {
    // 입출력 속도 향상
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    int n, m;
    cin >> n >> m;

    // 초기화: 자기 자신을 부모로 설정
    for (int i = 0; i <= n; i++) {
        parent[i] = i;
    }

    for (int i = 0; i < m; i++) {
        int op, a, b;
        cin >> op >> a >> b;

        if (op == 0) {
            // 합집합 연산
            unite(a, b);
        } else {
            // 같은 집합인지 확인 연산
            if (find(a) == find(b)) {
                cout << "YES\n";
            } else {
                cout << "NO\n";
            }
        }
    }

    return 0;
}