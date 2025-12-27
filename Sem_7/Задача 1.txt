#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cstdlib>

class Trajectory {
public:
    std::vector<double> t;
    std::vector<double> x;

    bool loadFromFile(const std::string& filename) {
        t.clear();
        x.clear();

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

            double tt = 0.0, xx = 0.0;
            if (!toDouble(a, tt)) continue;
            if (!toDouble(b, xx)) continue;

            t.push_back(tt);
            x.push_back(xx);
        }

        return t.size() >= 2 && x.size() == t.size();
    }

    void loadFromDefaults() {
        t.clear();
        x.clear();

        double t_arr[] = {0, 1, 2, 3, 4};
        double x_arr[] = {0, 2, 3.8, 6.1, 8.0};
        int n = 5;

        for (int i = 0; i < n; ++i) {
            t.push_back(t_arr[i]);
            x.push_back(x_arr[i]);
        }
    }

    std::vector<double> computeVelocity() const {
        std::vector<double> v;
        if (t.size() < 2) return v;

        v.reserve(t.size() - 1);
        for (size_t i = 0; i + 1 < t.size(); ++i) {
            double dt = t[i + 1] - t[i];
            if (dt == 0.0) v.push_back(0.0);
            else v.push_back((x[i + 1] - x[i]) / dt);
        }
        return v;
    }

    bool saveUsedCSV(const std::string& filename) const {
        std::ofstream out(filename.c_str());
        if (!out.is_open()) return false;

        out << "t,x\n";
        out << std::fixed << std::setprecision(6);
        for (size_t i = 0; i < t.size(); ++i) {
            out << t[i] << "," << x[i] << "\n";
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
        out << "set terminal pngcairo size 1000,600\n";
        out << "set output '" << out_png << "'\n";
        out << "set title 'x(t)'\n";
        out << "set xlabel 't'\n";
        out << "set ylabel 'x'\n";
        out << "plot '" << csv_name << "' using 1:2 with linespoints title 'x(t)'\n";
        return true;
    }

private:
    static std::string trim(const std::string& s) {
        size_t a = 0;
        while (a < s.size() && (s[a] == ' ' || s[a] == '\t' || s[a] == '\r' || s[a] == '\n')) a++;
        size_t b = s.size();
        while (b > a && (s[b - 1] == ' ' || s[b - 1] == ' ' || s[b - 1] == '\t' || s[b - 1] == '\r' || s[b - 1] == '\n')) b--;
        return s.substr(a, b - a);
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
};

int main() {
    Trajectory tr;

    bool ok = tr.loadFromFile("traj.csv");
    if (!ok) {
        tr.loadFromDefaults();
        std::cout << "traj.csv not found -> using default data in code\n";
    } else {
        std::cout << "traj.csv found -> using file data\n";
    }

    std::vector<double> v = tr.computeVelocity();

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Loaded points: " << tr.t.size() << "\n";
    std::cout << "Velocities:\n";
    for (size_t i = 0; i < v.size(); ++i) {
        std::cout << "v[" << i << "]=" << v[i] << "\n";
    }

    tr.saveUsedCSV("traj_used.csv");
    tr.generateGnuplotScript("plot.plt", "traj_used.csv", "plot.png");

    std::cout << "Saved: traj_used.csv, plot.plt\n";
    std::cout << "To build plot: gnuplot plot.plt\n";
    return 0;
}
