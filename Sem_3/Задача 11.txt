#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>

using namespace std;

int main() {
    double v = 20.0;
    double dt = 1.0;
    double T = 10.0;

    vector<double> times;
    vector<double> heights;

    double h = 0.0;

    for (double t = 0.0; t <= T; t += dt) {
        times.push_back(t);
        heights.push_back(h);
        h += v * dt;
    }

    vector<double>::iterator it =
        max_element(heights.begin(), heights.end());

    size_t idx = distance(heights.begin(), it);
    double h_max = heights[idx];
    double t_max = times[idx];

    cout << fixed << setprecision(2);
    cout << "Maximum height: " << h_max << endl;
    cout << "Time of maximum height: " << t_max << endl;

    return 0;
}
