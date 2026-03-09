#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    string line;
    vector<int> alphabet(26, 0);
    int max_freq = 0;

    while (getline(cin, line)) {
        for (char c : line) {
            if (c >= 'a' && c <= 'z') {
                int index = c - 'a';
                alphabet[index]++;
                if (alphabet[index] > max_freq) {
                    max_freq = alphabet[index];
                }
            }
        }
    }

    for (int i = 0; i < 26; ++i) {
        if (alphabet[i] == max_freq) {
            cout << (char)(i + 'a');
        }
    }
    cout << endl;

    return 0;
}