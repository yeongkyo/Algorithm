#include <iostream>
#include <string>
#include <cmath>
#include <iomanip>

using namespace std;

int days_in_month[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
string m_names[] = {"jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"};

int get_month_index(const string& m) {
    for (int i = 0; i < 12; ++i) {
        if (m_names[i] == m) return i;
    }
    return 0;
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    double l;
    int n;
    if (!(cin >> l >> n)) return 0;

    double pi = acos(-1.0);
    double epsilon = 23.439281 * pi / 180.0; 
    double phi = l * pi / 180.0;            

    int solstice_days = 0;
    for (int i = 0; i < 5; ++i) {
        solstice_days += days_in_month[i]; 
    }
    solstice_days += 20;
    int t_solstice = solstice_days * 24 + 12; 

    for (int i = 0; i < n; ++i) {
        int d, h;
        string m;
        cin >> d >> m >> h;

        int m_idx = get_month_index(m);
        int q_days = 0;
        for (int j = 0; j < m_idx; ++j) {
            q_days += days_in_month[j];
        }
        q_days += (d - 1);
        int t_q = q_days * 24 + h;

        int delta_t = t_q - t_solstice;

        double theta = (double)delta_t / 8760.0 * 2.0 * pi;
        double lambda = theta + (h - 12.0) * pi / 12.0;

        double sx = cos(epsilon) * cos(theta);
        double sy = sin(theta);
        double sz = sin(epsilon) * cos(theta);

        double ox = cos(phi) * cos(lambda);
        double oy = cos(phi) * sin(lambda);
        double oz = sin(phi);

        double dot_product = sx * ox + sy * oy + sz * oz;

        if (dot_product <= 0) {
            cout << fixed << setprecision(10) << 0.0 << "\n";
        } else {
            double alpha = asin(dot_product) * 180.0 / pi;
            cout << fixed << setprecision(10) << alpha << "\n";
        }
    }

    return 0;
}