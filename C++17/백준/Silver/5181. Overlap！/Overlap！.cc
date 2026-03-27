#include <bits/stdc++.h>
using namespace std;

struct Exam {
    string day;
    int start, finish;
};

int toMinute(const string& s) {
    return (s[0] - '0') * 600 + (s[1] - '0') * 60 + (s[3] - '0') * 10 + (s[4] - '0');
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int K;
    cin >> K;

    for (int tc = 1; tc <= K; tc++) {
        int m, n;
        cin >> m >> n;

        map<string, Exam> mp;

        for (int i = 0; i < m; i++) {
            string course, day, time;
            cin >> course >> day >> time;
            int dash = time.find('-');
            string s1 = time.substr(0, dash);
            string s2 = time.substr(dash + 1);
            mp[course] = {day, toMinute(s1), toMinute(s2)};
        }

        cin.ignore();

        int ans = 0;

        for (int i = 0; i < n; i++) {
            string line;
            getline(cin, line);
            stringstream ss(line);
            vector<string> courses;
            string c;
            while (ss >> c) courses.push_back(c);

            bool overlap = false;
            for (int a = 0; a < (int)courses.size() && !overlap; a++) {
                for (int b = a + 1; b < (int)courses.size(); b++) {
                    Exam e1 = mp[courses[a]];
                    Exam e2 = mp[courses[b]];
                    if (e1.day != e2.day) continue;
                    if (max(e1.start, e2.start) < min(e1.finish, e2.finish)) {
                        overlap = true;
                        break;
                    }
                }
            }

            if (overlap) ans++;
        }

        cout << "Data Set " << tc << ":\n";
        cout << ans << '\n';
    }

    return 0;
}