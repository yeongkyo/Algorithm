#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

void countingSort(vector<int>& v, int exp) {
    int n = v.size();
    vector<int> output(n);
    int count[10] = {0};

    for (int i = 0; i < n; i++) {
        count[(v[i] / exp) % 10]++;
    }

    for (int i = 1; i < 10; i++) {
        count[i] += count[i - 1];
    }

    for (int i = n - 1; i >= 0; i--) {
        output[count[(v[i] / exp) % 10] - 1] = v[i];
        count[(v[i] / exp) % 10]--;
    }

    for (int i = 0; i < n; i++) {
        v[i] = output[i];
    }
}

void radixSort(vector<int>& v) {
    int max_val = *max_element(v.begin(), v.end());

    for (int exp = 1; max_val / exp > 0; exp *= 10) {
        countingSort(v, exp);
    }
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int N;
    cin >> N;

    vector<int> v(N);
    for (int i = 0; i < N; i++) {
        cin >> v[i];
        v[i] += 1000000;
    }

    radixSort(v);

    for (int i = 0; i < N; i++) {
        cout << v[i] - 1000000 << "\n";
    }

    return 0;
}