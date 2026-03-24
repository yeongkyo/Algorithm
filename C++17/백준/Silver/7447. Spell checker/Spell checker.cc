#include <iostream>
#include <vector>
#include <string>
#include <cmath>

using namespace std;

bool is_similar(const string& a, const string& b) {
    int len_a = a.length();
    int len_b = b.length();
    if (abs(len_a - len_b) > 1) return false;

    if (len_a == len_b) {
        int diff = 0;
        for (int i = 0; i < len_a; ++i) {
            if (a[i] != b[i]) {
                diff++;
                if (diff > 1) return false;
            }
        }
        return diff == 1;
    }

    const string& shorter = len_a < len_b ? a : b;
    const string& longer = len_a < len_b ? b : a;

    int i = 0, j = 0, diff = 0;
    while (i < shorter.length() && j < longer.length()) {
        if (shorter[i] == longer[j]) {
            i++;
            j++;
        } else {
            diff++;
            if (diff > 1) return false;
            j++;
        }
    }
    return true;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    vector<string> dict;
    string word;
    while (cin >> word && word != "#") {
        dict.push_back(word);
    }

    while (cin >> word && word != "#") {
        bool correct = false;
        for (const string& d : dict) {
            if (d == word) {
                correct = true;
                break;
            }
        }

        if (correct) {
            cout << word << " is correct\n";
        } else {
            cout << word << ":";
            for (const string& d : dict) {
                if (is_similar(word, d)) {
                    cout << " " << d;
                }
            }
            cout << "\n";
        }
    }

    return 0;
}