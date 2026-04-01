#include <iostream>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int N;
    cin >> N;

    int count = N / 4;

    for (int i = 0; i < count; i++) {
        cout << "long ";
    }
    cout << "int\n";

    return 0;
}