#include <iostream>

using namespace std;

int main() {
    int n;
    cin >> n;

    int solved[12] = {0, 12, 11, 11, 10, 9, 9, 9, 8, 7, 6, 6};
    int penalty[12] = {0, 1600, 894, 1327, 1311, 1004, 1178, 1357, 837, 1055, 556, 773};

    cout << solved[n] << " " << penalty[n] << "\n";

    return 0;
}