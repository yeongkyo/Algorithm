#include <iostream>
#include <string>
#include <vector>
#include <cctype>

using namespace std;

string sanitize(const string& s) {
    string res = "";
    for (char c : s) {
        if (isalpha(c)) {
            res += tolower(c);
        }
    }
    return res;
}

vector<int> get_permutation(int len, const string& key) {
    vector<int> p(len, -1);
    int next_source = 0;
    int key_idx = 0;
    int dir = 1;

    while (next_source < len) {
        int step = (key_idx < key.length()) ? (key[key_idx] - 'a' + 2) : 1;
        int count = 0;

        if (dir == 1) {
            for (int i = 0; i < len; ++i) {
                if (p[i] == -1) {
                    count++;
                    if (count == step) {
                        p[i] = next_source++;
                        count = 0;
                        if (next_source == len) break;
                    }
                }
            }
        } else {
            for (int i = len - 1; i >= 0; --i) {
                if (p[i] == -1) {
                    count++;
                    if (count == step) {
                        p[i] = next_source++;
                        count = 0;
                        if (next_source == len) break;
                    }
                }
            }
        }
        dir *= -1;
        key_idx++;
    }
    return p;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    char t;
    string key;
    cin >> t >> key;
    cin.ignore();

    string text;
    getline(cin, text);

    if (t == 'E') {
        string plaintext = sanitize(text);
        int len = plaintext.length();
        vector<int> p = get_permutation(len, key);
        string ciphertext(len, ' ');
        for (int i = 0; i < len; ++i) {
            ciphertext[i] = plaintext[p[i]];
        }
        cout << ciphertext << endl;
    } else {
        string ciphertext = text;
        int len = ciphertext.length();
        vector<int> p = get_permutation(len, key);
        string plaintext(len, ' ');
        for (int i = 0; i < len; ++i) {
            plaintext[p[i]] = ciphertext[i];
        }
        cout << plaintext << endl;
    }

    return 0;
}