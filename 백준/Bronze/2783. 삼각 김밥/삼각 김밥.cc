#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>

using namespace std;

int main() {
    double x, y;
    cin >> x >> y;

    double min_price = (x / y) * 1000.0;

    int n;
    cin >> n;

    for (int i = 0; i < n; ++i) {
        double xi, yi;
        cin >> xi >> yi;
        double current_price = (xi / yi) * 1000.0;
        if (current_price < min_price) {
            min_price = current_price;
        }
    }

    cout << fixed << setprecision(2) << min_price << endl;

    return 0;
}