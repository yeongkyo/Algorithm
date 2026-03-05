#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <iomanip>

using namespace std;

int n;

double eval_Y(double X, const vector<double>& a, const vector<double>& b) {
    vector<double> y(n);
    for (int i = 0; i < n; ++i) {
        y[i] = a[i] * X + b[i];
    }
    
    nth_element(y.begin(), y.begin() + n / 2, y.end());
    double median = y[n / 2];
    
    double sum = 0;
    for (int i = 0; i < n; ++i) {
        sum += abs(y[i] - median);
    }
    return sum;
}

double ternary_search(const vector<double>& a, const vector<double>& b) {
    double l = -1.0, r = 1.0;
    for (int iter = 0; iter < 100; ++iter) {
        double m1 = l + (r - l) / 3.0;
        double m2 = r - (r - l) / 3.0;
        
        if (eval_Y(m1, a, b) < eval_Y(m2, a, b)) {
            r = m2;
        } else {
            l = m1;
        }
    }
    return eval_Y(l, a, b);
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    if (!(cin >> n)) return 0;

    vector<double> a(n), b(n);
    for (int i = 0; i < n; ++i) {
        cin >> a[i] >> b[i];
    }

    double ans1 = ternary_search(a, b);
    
    double ans2 = ternary_search(b, a);

    cout << fixed << setprecision(15) << min(ans1, ans2) << "\n";

    return 0;
}