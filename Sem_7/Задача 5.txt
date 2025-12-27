#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstdio>

class MotionAnalyzer {
public:
    std::vector<double> t;
    std::vector<double> x;
    std::vector<double> v;
    std::vector<double> a;

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
            std::string s1, s2;
            if (!std::getline(ss, s1, ',')) continue;
            if (!std::getline(ss, s2, ',')) continue;

            double tt = 0.0, xx = 0.0;
            if (!toDouble(s1, tt)) continue;
            if (!toDouble(s2, xx)) continue;

            t.push_back(tt);
            x.push_back(xx);
        }

        return t.size() >= 2 && x.size() == t.size();
    }

    void loadDefaults() {
        clear();
        double t_arr[] = {0,1,2,3,4,5,6};
        double x_arr[] = {0,2.1,4.5,7.2,10.4,14.1,18.5};
        int n = 7;

        for (int i = 0; i < n; ++i) {
            t.push_back(t_arr[i]);
            x.push_back(x_arr[i]);
        }
    }

    void computeVelocity() {
        v.assign(t.size(), 0.0);
        if (t.size() < 2) return;

        auto vel = [](double x0, double x1, double t0, double t1) {
            double dt = t1 - t0;
            if (dt == 0.0) return 0.0;
            return (x1 - x0) / dt;
        };

        for (size_t i = 0; i + 1 < t.size(); ++i) {
            v[i] = vel(x[i], x[i + 1], t[i], t[i + 1]);
        }

        v.back() = v.size() >= 2 ? v[v.size() - 2] : 0.0;
    }

    void computeAcceleration() {
        a.assign(t.size(), 0.0);
        if (t.size() < 2) return;

        auto acc = [](double v0, double v1, double t0, double t1) {
            double dt = t1 - t0;
            if (dt == 0.0) return 0.0;
            return (v1 - v0) / dt;
        };

        for (size_t i = 0; i + 1 < t.size(); ++i) {
            a[i] = acc(v[i], v[i + 1], t[i], t[i + 1]);
        }

        a.back() = a.size() >= 2 ? a[a.size() - 2] : 0.0;
    }

    bool saveResults(const std::string& filename) const {
        std::ofstream out(filename.c_str());
        if (!out.is_open()) return false;

        out << "t,x,v,a\n";
        out << std::fixed << std::setprecision(6);

        for (size_t i = 0; i < t.size(); ++i) {
            double vv = (i < v.size()) ? v[i] : 0.0;
            double aa = (i < a.size()) ? a[i] : 0.0;
            out << t[i] << "," << x[i] << "," << vv << "," << aa << "\n";
        }
        return true;
    }

    bool plotWithGnuplot(const std::string& csv, const std::string& out_png) const {
        FILE* gp = popen("gnuplot -persistent", "w");
        if (!gp) return false;

        std::fprintf(gp, "set datafile separator ','\n");
        std::fprintf(gp, "set grid\n");
        std::fprintf(gp, "set terminal pngcairo size 1100,650\n");
        std::fprintf(gp, "set output '%s'\n", out_png.c_str());
        std::fprintf(gp, "set title 'Velocity and Acceleration'\n");
        std::fprintf(gp, "set xlabel 't'\n");
        std::fprintf(gp, "plot '%s' using 1:3 with linespoints title 'v(t)', \\\n", csv.c_str());
        std::fprintf(gp, "     '%s' using 1:4 with linespoints title 'a(t)'\n", csv.c_str());

        int rc = pclose(gp);
        return rc == 0;
    }

private:
    void clear() {
        t.clear();
        x.clear();
        v.clear();
        a.clear();
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
    MotionAnalyzer ma;

    bool ok = ma.loadFromFile("motion.csv");
    if (!ok) {
        ma.loadDefaults();
        std::cout << "motion.csv not found -> using default data in code\n";
    } else {
        std::cout << "motion.csv found -> using file data\n";
    }

    ma.computeVelocity();
    ma.computeAcceleration();

    if (!ma.saveResults("motion_processed.csv")) {
        std::cout << "Error: cannot write motion_processed.csv\n";
        return 1;
    }

    std::cout << "Saved: motion_processed.csv\n";

    bool plotted = ma.plotWithGnuplot("motion_processed.csv", "plot.png");
    if (plotted) std::cout << "Plot saved: plot.png\n";
    else std::cout << "Gnuplot not executed (check gnuplot installation)\n";

    return 0;
}
