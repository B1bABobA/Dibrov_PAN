#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>

class TelemetryFilter {
private:
    std::vector<std::vector<std::string>> data;
    int removed_ = 0;
    int kept_ = 0;

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

    static std::vector<std::string> splitCSV(const std::string& line) {
        std::vector<std::string> parts;
        std::stringstream ss(line);
        std::string token;
        while (std::getline(ss, token, ',')) {
            parts.push_back(trim(token));
        }
        return parts;
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

    static bool isValidAltitude(double alt) {
        return alt >= 0.0 && alt <= 20000.0;
    }

    static bool isValidSpeed(double spd) {
        return spd >= 0.0 && spd <= 500.0;
    }

public:
    bool loadFromCSV(const std::string& filename) {
        data.clear();
        removed_ = 0;
        kept_ = 0;

        std::ifstream in(filename.c_str());
        if (!in.is_open()) return false;

        std::string line;
        while (std::getline(in, line)) {
            line = stripComment(line);
            line = trim(line);
            if (line.empty()) continue;

            std::vector<std::string> row = splitCSV(line);
            if (!row.empty()) data.push_back(row);
        }
        return !data.empty();
    }

    void filterData() {
        if (data.empty()) return;

        std::vector<std::vector<std::string>> out;
        out.push_back(data[0]);

        for (size_t i = 1; i < data.size(); ++i) {
            const std::vector<std::string>& row = data[i];
            if (row.size() < 5) {
                removed_++;
                continue;
            }

            double alt = 0.0, spd = 0.0;
            bool okAlt = toDouble(row[1], alt);
            bool okSpd = toDouble(row[2], spd);

            if (!okAlt || !okSpd) {
                removed_++;
                continue;
            }

            if (!isValidAltitude(alt) || !isValidSpeed(spd)) {
                removed_++;
                continue;
            }

            out.push_back(row);
            kept_++;
        }

        data = out;
    }

    bool saveToCSV(const std::string& filename) {
        std::ofstream out(filename.c_str());
        if (!out.is_open()) return false;

        for (size_t i = 0; i < data.size(); ++i) {
            for (size_t j = 0; j < data[i].size(); ++j) {
                if (j) out << ",";
                out << data[i][j];
            }
            out << "\n";
        }
        return true;
    }

    void printFilteredStats() {
        int total_rows = 0;
        if (!data.empty()) {
            total_rows = (int)data.size() - 1;
        }
        std::cout << "Telemetry filter stats\n";
        std::cout << "Kept rows: " << kept_ << "\n";
        std::cout << "Removed rows: " << removed_ << "\n";
        std::cout << "Rows in output (without header): " << total_rows << "\n";
    }
};

int main() {
    TelemetryFilter tf;

    if (!tf.loadFromCSV("telemetry.csv")) {
        std::cout << "Error: cannot open telemetry.csv\n";
        return 1;
    }

    tf.filterData();

    if (!tf.saveToCSV("telemetry_filtered.csv")) {
        std::cout << "Error: cannot write telemetry_filtered.csv\n";
        return 1;
    }

    tf.printFilteredStats();
    std::cout << "Saved to telemetry_filtered.csv\n";
    return 0;
}
