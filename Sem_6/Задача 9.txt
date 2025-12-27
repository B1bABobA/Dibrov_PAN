#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <cstdlib>

class Aircraft {
public:
    bool loadFromFile(const std::string& filename) {
        std::ifstream in(filename.c_str());
        if (!in.is_open()) return false;

        std::string line;
        while (std::getline(in, line)) {
            line = trim(stripComment(line));
            if (line.empty()) continue;

            size_t eq = line.find('=');
            if (eq == std::string::npos) continue;

            std::string key = trim(line.substr(0, eq));
            std::string val = trim(line.substr(eq + 1));

            double v = 0.0;
            if (!toDouble(val, v)) continue;

            if (key == "mass") mass = v;
            else if (key == "wing_area") wing_area = v;
            else if (key == "Cx") Cx = v;
            else if (key == "Cy") Cy = v;
            else if (key == "max_thrust") max_thrust = v;
            else if (key == "specific_consumption") specific_consumption = v;
        }

        return (mass > 0 && wing_area > 0 && max_thrust > 0);
    }

    double mass = 0.0;
    double wing_area = 0.0;
    double Cx = 0.0;
    double Cy = 0.0;
    double max_thrust = 0.0;
    double specific_consumption = 0.0;

private:
    static std::string stripComment(const std::string& s) {
        size_t p = s.find("//");
        if (p == std::string::npos) return s;
        return s.substr(0, p);
    }

    static std::string trim(const std::string& s) {
        size_t a = 0;
        while (a < s.size() && (s[a] == ' ' || s[a] == '\t' || s[a] == '\r' || s[a] == '\n')) a++;
        size_t b = s.size();
        while (b > a && (s[b - 1] == ' ' || s[b - 1] == '\t' || s[b - 1] == '\r' || s[b - 1] == '\n')) b--;
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

class Environment {
public:
    bool loadAtmosphereTable(const std::string& filename) {
        alt.clear();
        dens.clear();

        std::ifstream in(filename.c_str());
        if (!in.is_open()) return false;

        std::string line;
        if (!std::getline(in, line)) return false;

        while (std::getline(in, line)) {
            line = trim(stripComment(line));
            if (line.empty()) continue;

            std::vector<std::string> parts = splitCSV(line);
            if (parts.size() < 2) continue;

            double a = 0.0, d = 0.0;
            if (!toDouble(parts[0], a)) continue;
            if (!toDouble(parts[1], d)) continue;

            alt.push_back(a);
            dens.push_back(d);
        }

        return alt.size() >= 2;
    }

    double getDensity(double altitude) {
        if (alt.empty()) return 0.0;
        if (altitude <= alt.front()) return dens.front();
        if (altitude >= alt.back()) return dens.back();

        for (size_t i = 1; i < alt.size(); ++i) {
            if (altitude <= alt[i]) {
                double a0 = alt[i - 1], a1 = alt[i];
                double d0 = dens[i - 1], d1 = dens[i];
                double t = (altitude - a0) / (a1 - a0);
                return d0 + (d1 - d0) * t;
            }
        }
        return dens.back();
    }

private:
    std::vector<double> alt;
    std::vector<double> dens;

    static std::string stripComment(const std::string& s) {
        size_t p = s.find("//");
        if (p == std::string::npos) return s;
        return s.substr(0, p);
    }

    static std::string trim(const std::string& s) {
        size_t a = 0;
        while (a < s.size() && (s[a] == ' ' || s[a] == '\t' || s[a] == '\r' || s[a] == '\n')) a++;
        size_t b = s.size();
        while (b > a && (s[b - 1] == ' ' || s[b - 1] == '\t' || s[b - 1] == '\r' || s[b - 1] == '\n')) b--;
        return s.substr(a, b - a);
    }

    static std::vector<std::string> splitCSV(const std::string& line) {
        std::vector<std::string> parts;
        std::stringstream ss(line);
        std::string token;
        while (std::getline(ss, token, ',')) parts.push_back(trim(token));
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
};

int main() {
    Aircraft ac;
    Environment env;

    if (!ac.loadFromFile("aircraft_params.txt")) {
        std::cout << "Error: cannot load aircraft_params.txt\n";
        return 1;
    }
    if (!env.loadAtmosphereTable("atmosphere.csv")) {
        std::cout << "Error: cannot load atmosphere.csv\n";
        return 1;
    }

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Aircraft loaded\n";
    std::cout << "mass=" << ac.mass << "\n";
    std::cout << "wing_area=" << ac.wing_area << "\n";
    std::cout << "Cx=" << ac.Cx << "\n";
    std::cout << "Cy=" << ac.Cy << "\n";
    std::cout << "max_thrust=" << ac.max_thrust << "\n";
    std::cout << "specific_consumption=" << ac.specific_consumption << "\n\n";

    double h1 = 500.0;
    double h2 = 1500.0;

    std::cout << "Density interpolation\n";
    std::cout << "Altitude " << h1 << " m: density=" << env.getDensity(h1) << "\n";
    std::cout << "Altitude " << h2 << " m: density=" << env.getDensity(h2) << "\n";

    return 0;
}
