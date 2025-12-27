#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <utility>
#include <algorithm>
#include <iomanip>
#include <cstdlib>

class AltitudeFilter {
public:
    std::vector<std::pair<double,double>> data;
    std::vector<std::pair<double,double>> filtered;

    bool loadFromFile(const std::string& filename) {
        data.clear();
        filtered.clear();

        std::ifstream in(filename.c_str());
        if (!in.is_open()) return false;

        std::string line;
        if (!std::getline(in, line)) return false;

        while (std::getline(in, line)) {
            line = trim(line);
            if (line.empty()) continue;

            std::stringstream ss(line);
            std::string a, b;
            if (!std::getline(ss, a, ',')) continue;
            if (!std::getline(ss, b, ',')) continue;

            double t = 0.0, h = 0.0;
            if (!toDouble(a, t)) continue;
            if (!toDouble(b, h)) continue;

            data.push_back(std::make_pair(t, h));
        }

        return data.size() >= 1;
    }

    void loadDefaults() {
        data.clear();
        filtered.clear();

        double t_arr[] = {0,1,2,3,4,5};
        double h_arr[] = {1000,1005,1010,1500,1020,1025};
        int n = 6;

        for (int i = 0; i < n; ++i) {
            data.push_back(std::make_pair(t_arr[i], h_arr[i]));
        }
    }

    void filterWithLambda() {
        filtered = data;

        filtered.erase(
            std::remove_if(filtered.begin(), filtered.end(),
                [](const std::pair<double,double>& p) {
                    double H = p.second;
                    return !(H > 900.0 && H < 1100.0);
                }
            ),
            filtered.end()
        );
    }

    bool saveCSV(const std::string& filename, const std::vector<std::pair<double,double>>& v) const {
        std::ofstream out(filename.c_str());
        if (!out.is_open()) return false;

        out << "t,H\n";
        out << std::fixed << std::setprecision(6);
        for (size_t i = 0; i < v.size(); ++i) {
            out << v[i].first << "," << v[i].second << "\n";
        }
        return true;
    }

    bool generateGnuplotScript(const std::string& script_name,
                              const std::string& orig_csv,
                              const std::string& filt_csv,
                              const std::string& out_png) const {
        std::ofstream out(script_name.c_str());
        if (!out.is_open()) return false;

        out << "set datafile separator ','\n";
        out << "set grid\n";
        out << "set terminal pngcairo size 1100,650\n";
        out << "set output '" << out_png << "'\n";
        out << "set title 'Altitude: original vs filtered'\n";
        out << "set xlabel 't'\n";
        out << "set ylabel 'H'\n";
        out << "plot '" << orig_csv << "' using 1:2 with linespoints title 'original',\\\n";
        out << "     '" << filt_csv << "' using 1:2 with linespoints title 'filtered'\n";
        return true;
    }

private:
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
    AltitudeFilter af;

    bool ok = af.loadFromFile("altitude.csv");
    if (!ok) {
        af.loadDefaults();
        std::cout << "altitude.csv not found -> using default data in code\n";
    } else {
        std::cout << "altitude.csv found -> using file data\n";
    }

    af.filterWithLambda();

    if (!af.saveCSV("original.csv", af.data)) {
        std::cout << "Error: cannot write original.csv\n";
        return 1;
    }
    if (!af.saveCSV("filtered.csv", af.filtered)) {
        std::cout << "Error: cannot write filtered.csv\n";
        return 1;
    }
    if (!af.generateGnuplotScript("plot.plt", "original.csv", "filtered.csv", "plot.png")) {
        std::cout << "Error: cannot write plot.plt\n";
        return 1;
    }

    std::cout << "Saved: original.csv, filtered.csv, plot.plt\n";
    std::cout << "To build plot: gnuplot plot.plt\n";
    return 0;
}
