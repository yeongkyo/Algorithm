#include <iostream>
#include <vector>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int A, B;
    if (cin >> A >> B) {
        int m;
        cin >> m;

        long long decimalValue = 0;
        for (int i = 0; i < m; ++i) {
            int digit;
            cin >> digit;
            decimalValue = decimalValue * A + digit;
        }

        if (decimalValue == 0) {
            cout << 0;
        } else {
            vector<int> result;
            while (decimalValue > 0) {
                result.push_back(decimalValue % B);
                decimalValue /= B;
            }

            for (int i = result.size() - 1; i >= 0; --i) {
                cout << result[i] << (i > 0 ? " " : "");
            }
        }
    }
    return 0;
}