#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>
#include <cstdlib>

class Navigator {
public:
    std::vector<double> t;
    std::vector<double> x;
    std::vector<double> y;

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

            double tt = 0.0, xx = 0.0, yy = 0.0;
            if (!toDouble(a, tt)) continue;
            if (!toDouble(b, xx)) continue;
            if (!toDouble(c, yy)) continue;

            t.push_back(tt);
            x.push_back(xx);
            y.push_back(yy);
        }

        return t.size() >= 2 && x.size() == t.size() && y.size() == t.size();
    }

    void loadDefaults() {
        clear();
        double t_arr[] = {0,1,2,3,4};
        double x_arr[] = {0,5,9,12,14};
        double y_arr[] = {0,1,3,6,10};
        int n = 5;

        for (int i = 0; i < n; ++i) {
            t.push_back(t_arr[i]);
            x.push_back(x_arr[i]);
            y.push_back(y_arr[i]);
        }
    }

    std::vector<double> computeSpeedMagnitude() const {
        std::vector<double> v;
        if (t.size() < 2) return v;

        v.reserve(t.size() - 1);
        for (size_t i = 0; i + 1 < t.size(); ++i) {
            double dt = t[i + 1] - t[i];
            if (dt == 0.0) {
                v.push_back(0.0);
            } else {
                double vx = (x[i + 1] - x[i]) / dt;
                double vy = (y[i + 1] - y[i]) / dt;
                v.push_back(std::sqrt(vx * vx + vy * vy));
            }
        }
        return v;
    }

    std::vector<double> timeForSpeed() const {
        std::vector<double> tt;
        if (t.size() < 2) return tt;
        tt.reserve(t.size() - 1);
        for (size_t i = 0; i + 1 < t.size(); ++i) {
            tt.push_back(t[i]);
        }
        return tt;
    }

    bool saveTV(const std::string& filename, const std::vector<double>& tt, const std::vector<double>& v) const {
        if (tt.size() != v.size()) return false;

        std::ofstream out(filename.c_str());
        if (!out.is_open()) return false;

        out << "t,v\n";
        out << std::fixed << std::setprecision(6);
        for (size_t i = 0; i < tt.size(); ++i) {
            out << tt[i] << "," << v[i] << "\n";
        }
        return true;
    }

    bool generatePlot(const std::string& script_name, const std::string& tv_csv, const std::string& out_png) const {
        std::ofstream out(script_name.c_str());
        if (!out.is_open()) return false;

        out << "set datafile separator ','\n";
        out << "set grid\n";
        out << "set terminal pngcairo size 1100,650\n";
        out << "set output '" << out_png << "'\n";
        out << "set title 'Speed magnitude v(t)'\n";
        out << "set xlabel 't'\n";
        out << "set ylabel 'v'\n";
        out << "plot '" << tv_csv << "' using 1:2 with linespoints title 'v(t)'\n";
        return true;
    }

private:
    void clear() {
        t.clear();
        x.clear();
        y.clear();
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
    Navigator nav;

    bool ok = nav.loadFromFile("navigation.csv");
    if (!ok) {
        nav.loadDefaults();
        std::cout << "navigation.csv not found -> using default data in code\n";
    } else {
        std::cout << "navigation.csv found -> using file data\n";
    }

    std::vector<double> tt = nav.timeForSpeed();
    std::vector<double> v = nav.computeSpeedMagnitude();

    if (!nav.saveTV("tv.csv", tt, v)) {
        std::cout << "Error: cannot write tv.csv\n";
        return 1;
    }

    if (!nav.generatePlot("plot.plt", "tv.csv", "plot.png")) {
        std::cout << "Error: cannot write plot.plt\n";
        return 1;
    }

    std::cout << "Saved: tv.csv, plot.plt\n";
    std::cout << "To build plot: gnuplot plot.plt\n";
    return 0;
}
