#include <iostream>
#include <cmath>

using namespace std;

int main() {
    int a, b;
    if (!(cin >> a >> b)) return 0;

    int x1 = (a - 1) / 4;
    int y1 = (a - 1) % 4;
    
    int x2 = (b - 1) / 4;
    int y2 = (b - 1) % 4;

    int distance = abs(x1 - x2) + abs(y1 - y2);

    cout << distance << endl;

    return 0;
}