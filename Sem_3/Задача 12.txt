#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

int main() {
    vector<double> sensorData;

    sensorData.push_back(0.0);
    sensorData.push_back(1.5);
    sensorData.push_back(0.0);
    sensorData.push_back(-0.3);
    sensorData.push_back(0.0);
    sensorData.push_back(2.1);

    int zeroCount = count(sensorData.begin(), sensorData.end(), 0.0);

    cout << "Number of zero measurements: " << zeroCount << endl;

    return 0;
}
