#include <iostream>
#include <string>
#include <vector>

using namespace std;

int get_holes(char c) {
    if (c == 'B') return 2;
    if (c == 'A' || c == 'D' || c == 'O' || c == 'P' || c == 'Q' || c == 'R') return 1;
    return 0;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    string C1, C2;
    if (!(cin >> C1 >> C2)) return 0;

    vector<int> holes;
    for (char c : C1) {
        holes.push_back(get_holes(c));
    }

    vector<string> lines(5, "");
    bool placed_c2 = false;

    for (size_t i = 0; i < holes.size(); ++i) {
        int h = holes[i];
        if (h == 0 && !placed_c2) {
            lines[0] += C2;
            string dots(C2.length(), '.');
            for (int j = 1; j < 5; ++j) {
                lines[j] += dots;
            }
            placed_c2 = true;
        } else if (h == 0) {
            for (int j = 0; j < 5; ++j) lines[j] += "Z";
        } else if (h == 1) {
            lines[0] += "ZZZ";
            lines[1] += "Z.Z";
            lines[2] += "ZZZ";
            lines[3] += "Z..";
            lines[4] += "Z..";
        } else if (h == 2) {
            lines[0] += "ZZZ";
            lines[1] += "Z.Z";
            lines[2] += "ZZZ";
            lines[3] += "Z.Z";
            lines[4] += "ZZZ";
        }

        if (i + 1 < holes.size()) {
            for (int j = 0; j < 5; ++j) lines[j] += ".";
        }
    }

    for (int i = 0; i < 5; ++i) {
        cout << lines[i] << "\n";
    }

    return 0;
}