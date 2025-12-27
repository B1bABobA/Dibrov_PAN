#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>

using namespace std;

int main() {
    vector<double> y;

    y.push_back(0.0);
    y.push_back(5.0);
    y.push_back(3.0);
    y.push_back(10.0);
    y.push_back(7.0);

    if (y.size() < 2) {
        cout << "Total path: 0" << endl;
        return 0;
    }

    vector<double> diff(y.size());

    adjacent_difference(y.begin(), y.end(), diff.begin());

    double totalPath = 0.0;

    for (size_t i = 1; i < diff.size(); i++) {
        totalPath += fabs(diff[i]);
    }

    cout << "Total path: " << totalPath << endl;

    return 0;
}
