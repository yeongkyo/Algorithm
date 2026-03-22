#include <iostream>
#include <string>

using namespace std;

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    string input;
    while (cin >> input && input != "end") {
        if (input == "animal") {
            cout << "Panthera tigris\n";
        } else if (input == "tree") {
            cout << "Pinus densiflora\n";
        } else if (input == "flower") {
            cout << "Forsythia koreana\n";
        }
    }

    return 0;
}