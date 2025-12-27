#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

int main() {
    vector<double> x;

    x.push_back(5.0);
    x.push_back(2.0);
    x.push_back(8.0);
    x.push_back(5.0);
    x.push_back(3.0);
    x.push_back(2.0);
    x.push_back(10.0);

    sort(x.begin(), x.end());

    vector<double>::iterator new_end = unique(x.begin(), x.end());
    x.erase(new_end, x.end());

    for (size_t i = 0; i < x.size(); i++) {
        cout << x[i] << " ";
    }
    cout << endl;

    return 0;
}
