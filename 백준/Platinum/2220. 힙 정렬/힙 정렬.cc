#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    if (!(cin >> n)) return 0;

    vector<int> heap(n + 1);
    
    heap[1] = 1;

    for (int i = 2; i <= n; ++i) {
        int idx = i - 1;
        
        while (idx > 1) {
            swap(heap[idx], heap[idx / 2]);
            idx /= 2;
        }
        
        heap[i] = heap[1];
        heap[1] = i;
    }

    for (int i = 1; i <= n; ++i) {
        cout << heap[i] << (i == n ? "" : " ");
    }
    cout << "\n";

    return 0;
}