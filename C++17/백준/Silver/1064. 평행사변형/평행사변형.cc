#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

using namespace std;

double get_dist(double x1, double y1, double x2, double y2) {
    return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2));
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);

    double xA, yA, xB, yB, xC, yC;
    if (!(cin >> xA >> yA >> xB >> yB >> xC >> yC)) return 0;

    double cross_product = (xB - xA) * (yC - yA) - (xC - xA) * (yB - yA);
    if (cross_product == 0) {
        cout << "-1.0" << endl;
        return 0;
    }

    double dAB = get_dist(xA, yA, xB, yB);
    double dBC = get_dist(xB, yB, xC, yC);
    double dCA = get_dist(xC, yC, xA, yA);

    vector<double> perimeters;
    perimeters.push_back(2 * (dAB + dBC));
    perimeters.push_back(2 * (dBC + dCA));
    perimeters.push_back(2 * (dCA + dAB));

    double max_p = *max_element(perimeters.begin(), perimeters.end());
    double min_p = *min_element(perimeters.begin(), perimeters.end());

    cout << fixed << setprecision(16) << max_p - min_p << endl;

    return 0;
}