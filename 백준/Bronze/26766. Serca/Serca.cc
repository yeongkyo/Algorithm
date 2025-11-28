#include <iostream>
#include <string>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    int n;
    cin >> n;

    string heart = 
        " @@@   @@@ \n"
        "@   @ @   @\n"
        "@    @    @\n"
        "@         @\n"
        " @       @ \n"
        "  @     @  \n"
        "   @   @   \n"
        "    @ @    \n"
        "     @     \n";

    for (int i = 0; i < n; ++i) {
        cout << heart;
    }

    return 0;
}