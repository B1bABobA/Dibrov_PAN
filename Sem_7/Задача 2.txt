#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <cstdlib>

class SensorData {
public:
    std::vector<double> t;
    std::vector<double> h1;
    std::vector<double> h2;
    std::vector<double> dh;

    bool loadFromFile(const std::string& filename) {
        clear();

        std::ifstream in(filename.c_str());
        if (!in.is_open()) return false;

        std::string line;
        if (!std::getline(in, line)) return false;

        while (std::getline(in, line)) {
            line = trim(line);
            if (line.empty()) continue;

            std::stringstream ss(line);
            std::string a, b, c;
            if (!std::getline(ss, a, ',')) continue;
            if (!std::getline(ss, b, ',')) continue;
            if (!std::getline(ss, c, ',')) continue;

            double tt = 0.0, v1 = 0.0, v2 = 0.0;
            if (!toDouble(a, tt)) continue;
            if (!toDouble(b, v1)) continue;
            if (!toDouble(c, v2)) continue;

            t.push_back(tt);
            h1.push_back(v1);
            h2.push_back(v2);
        }

        return t.size() >= 2 && h1.size() == t.size() && h2.size() == t.size();
    }

    void loadDefaults() {
        clear();

        double t_arr[]  = {0, 1, 2, 3, 4};
        double h1_arr[] = {1000, 1015, 1030, 1048, 1060};
        double h2_arr[] = {1002, 1012, 1028, 1042, 1058};
        int n = 5;

        for (int i = 0; i < n; ++i) {
            t.push_back(t_arr[i]);
            h1.push_back(h1_arr[i]);
            h2.push_back(h2_arr[i]);
        }
    }

    void computeDiff() {
        dh.clear();
        dh.reserve(t.size());

        auto diffFunc = [](double a, double b) {
            return std::fabs(a - b);
        };

        for (size_t i = 0; i < t.size(); ++i) {
            dh.push_back(diffFunc(h1[i], h2[i]));
        }
    }

    bool saveDiffCSV(const std::string& filename) const {
        std::ofstream out(filename.c_str());
        if (!out.is_open()) return false;

        out << "t,h1,h2,dh\n";
        out << std::fixed << std::setprecision(6);

        for (size_t i = 0; i < t.size(); ++i) {
            out << t[i] << "," << h1[i] << "," << h2[i] << "," << dh[i] << "\n";
        }
        return true;
    }

    bool generateGnuplotScript(const std::string& script_name,
                              const std::string& csv_name,
                              const std::string& out_png) const {
        std::ofstream out(script_name.c_str());
        if (!out.is_open()) return false;

        out << "set datafile separator ','\n";
        out << "set grid\n";
        out << "set terminal pngcairo size 1100,650\n";
        out << "set output '" << out_png << "'\n";
        out << "set title 'Sensors comparison'\n";
        out << "set xlabel 't'\n";
        out << "set ylabel 'h'\n";
        out << "plot '" << csv_name << "' using 1:2 with linespoints title 'h1(t)',\\\n";
        out << "     '" << csv_name << "' using 1:3 with linespoints title 'h2(t)',\\\n";
        out << "     '" << csv_name << "' using 1:4 with linespoints title 'dh(t)'\n";
        return true;
    }

private:
    void clear() {
        t.clear();
        h1.clear();
        h2.clear();
        dh.clear();
    }

    static std::string trim(const std::string& s) {
        size_t a = 0;
        while (a < s.size() && (s[a] == ' ' || s[a] == '\t' || s[a] == '\r' || s[a] == '\n')) a++;
        size_t b = s.size();
        while (b > a && (s[b - 1] == ' ' || s[b - 1] == '\t' || s[b - 1] == '\r' || s[b - 1] == '\n')) b--;
        return s.substr(a, b - a);
    }

    static bool toDouble(const std::string& s, double& out) {
        std::string tt = trim(s);
        if (tt.empty()) return false;
        char* endp = 0;
        out = std::strtod(tt.c_str(), &endp);
        if (endp == tt.c_str()) return false;
        while (*endp == ' ' || *endp == '\t' || *endp == '\r' || *endp == '\n') endp++;
        return *endp == '\0';
    }
};

int main() {
    SensorData sd;

    bool ok = sd.loadFromFile("sensors.csv");
    if (!ok) {
        sd.loadDefaults();
        std::cout << "sensors.csv not found -> using default data in code\n";
    } else {
        std::cout << "sensors.csv found -> using file data\n";
    }

    sd.computeDiff();

    if (!sd.saveDiffCSV("diff.csv")) {
        std::cout << "Error: cannot write diff.csv\n";
        return 1;
    }

    if (!sd.generateGnuplotScript("plot.plt", "diff.csv", "plot.png")) {
        std::cout << "Error: cannot write plot.plt\n";
        return 1;
    }

    std::cout << "Saved: diff.csv, plot.plt\n";
    std::cout << "To build plot: gnuplot plot.plt\n";
    return 0;
}
