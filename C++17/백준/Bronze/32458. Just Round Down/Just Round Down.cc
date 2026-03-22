#include <iostream>
#include <string>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    string s;
    cin >> s;

    size_t pos = s.find('.');
    cout << s.substr(0, pos) << "\n";

    return 0;
}