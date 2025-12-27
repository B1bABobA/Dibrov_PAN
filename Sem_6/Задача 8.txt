#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cstdlib>

class FuelAnalyzer {
private:
    std::vector<double> time_data;
    std::vector<double> fuel_data;
    std::vector<double> rpm_data;

    std::vector<int> anomaly_idx;

    static std::string trim(const std::string& s) {
        size_t a = 0;
        while (a < s.size() && (s[a] == ' ' || s[a] == '\t' || s[a] == '\r' || s[a] == '\n')) a++;
        size_t b = s.size();
        while (b > a && (s[b - 1] == ' ' || s[b - 1] == '\t' || s[b - 1] == '\r' || s[b - 1] == '\n')) b--;
        return s.substr(a, b - a);
    }

    static std::string stripComment(const std::string& s) {
        size_t p = s.find("//");
        if (p == std::string::npos) return s;
        return s.substr(0, p);
    }

    static bool toDouble(const std::string& s, double& out) {
        std::string t = trim(s);
        if (t.empty()) return false;
        char* endp = 0;
        out = std::strtod(t.c_str(), &endp);
        if (endp == t.c_str()) return false;
        while (*endp == ' ' || *endp == '\t' || *endp == '\r' || *endp == '\n') endp++;
        return *endp == '\0';
    }

    bool isAnomaly(double consumption, double avg, double threshold) {
        return consumption > avg * threshold;
    }

public:
    bool loadData(const std::string& filename) {
        time_data.clear();
        fuel_data.clear();
        rpm_data.clear();
        anomaly_idx.clear();

        std::ifstream in(filename.c_str());
        if (!in.is_open()) return false;

        std::string line;
        if (!std::getline(in, line)) return false;

        while (std::getline(in, line)) {
            line = stripComment(line);
            line = trim(line);
            if (line.empty()) continue;

            std::stringstream ss(line);
            std::string token;
            std::vector<std::string> parts;

            while (std::getline(ss, token, ',')) {
                parts.push_back(trim(token));
            }
            if (parts.size() != 3) continue;

            double t = 0.0, f = 0.0, r = 0.0;
            if (!toDouble(parts[0], t)) continue;
            if (!toDouble(parts[1], f)) continue;
            if (!toDouble(parts[2], r)) continue;

            time_data.push_back(t);
            fuel_data.push_back(f);
            rpm_data.push_back(r);
        }

        return !time_data.empty();
    }

    double calculateAverageConsumption() {
        if (fuel_data.empty()) return 0.0;
        double sum = 0.0;
        for (size_t i = 0; i < fuel_data.size(); ++i) sum += fuel_data[i];
        return sum / (double)fuel_data.size();
    }

    void detectAnomalies() {
        anomaly_idx.clear();
        double avg = calculateAverageConsumption();
        double threshold = 1.5;

        for (size_t i = 0; i < fuel_data.size(); ++i) {
            if (isAnomaly(fuel_data[i], avg, threshold)) {
                anomaly_idx.push_back((int)i);
            }
        }
    }

    void generateReport(const std::string& filename) {
        std::ofstream out(filename.c_str());
        if (!out.is_open()) return;

        double avg = calculateAverageConsumption();

        out << std::fixed << std::setprecision(1);
        out << "Fuel consumption report\n";
        out << "Records: " << fuel_data.size() << "\n";
        out << "Average consumption: " << avg << "\n";
        out << "Anomaly rule: consumption > avg * 1.5\n\n";

        if (anomaly_idx.empty()) {
            out << "No anomalies detected\n";
            return;
        }

        out << "Anomalies detected:\n";
        for (size_t k = 0; k < anomaly_idx.size(); ++k) {
            int i = anomaly_idx[k];
            out << "Time: " << time_data[i]
                << ", Consumption: " << fuel_data[i]
                << ", RPM: " << rpm_data[i] << "\n";
        }
    }

    void printSummary() {
        double avg = calculateAverageConsumption();
        std::cout << std::fixed << std::setprecision(1);
        std::cout << "Fuel analyzer summary\n";
        std::cout << "Records: " << fuel_data.size() << "\n";
        std::cout << "Average consumption: " << avg << "\n";
        std::cout << "Anomalies: " << anomaly_idx.size() << "\n";
        for (size_t k = 0; k < anomaly_idx.size(); ++k) {
            int i = anomaly_idx[k];
            std::cout << "Anomaly at time " << time_data[i]
                      << ": consumption=" << fuel_data[i]
                      << ", rpm=" << rpm_data[i] << "\n";
        }
    }
};

int main() {
    FuelAnalyzer fa;

    if (!fa.loadData("fuel_data.csv")) {
        std::cout << "Error: cannot open fuel_data.csv\n";
        return 1;
    }

    fa.detectAnomalies();
    fa.generateReport("fuel_report.txt");
    fa.printSummary();

    std::cout << "Saved to fuel_report.txt\n";
    return 0;
}
