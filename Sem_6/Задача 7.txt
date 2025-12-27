#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>
#include <iomanip>
#include <algorithm>

class WaypointSorter {
private:
    struct Waypoint {
        int id;
        double x, y, z;
        std::string name;
        double distance;
    };

    std::vector<Waypoint> waypoints;

    static std::string trim(const std::string& s) {
        size_t a = 0;
        while (a < s.size() && (s[a] == ' ' || s[a] == '\t' || s[a] == '\r' || s[a] == '\n')) a++;
        size_t b = s.size();
        while (b > a && (s[b - 1] == ' ' || s[b - 1] == '\t' || s[b - 1] == '\r' || s[b - 1] == '\n')) b--;
        return s.substr(a, b - a);
    }

    static bool compareByDistance(const Waypoint& a, const Waypoint& b) {
        return a.distance < b.distance;
    }

public:
    bool loadWaypoints(const std::string& filename) {
        waypoints.clear();

        std::ifstream in(filename.c_str());
        if (!in.is_open()) return false;

        std::string line;
        while (std::getline(in, line)) {
            line = trim(line);
            if (line.empty()) continue;

            std::stringstream ss(line);
            std::string token;
            std::vector<std::string> parts;

            while (std::getline(ss, token, ',')) {
                parts.push_back(trim(token));
            }

            if (parts.size() != 5) continue;

            Waypoint w;
            w.id = std::atoi(parts[0].c_str());
            w.x = std::atof(parts[1].c_str());
            w.y = std::atof(parts[2].c_str());
            w.z = std::atof(parts[3].c_str());
            w.name = parts[4];
            w.distance = 0.0;

            waypoints.push_back(w);
        }

        return !waypoints.empty();
    }

    void calculateDistances(double current_x, double current_y, double current_z) {
        for (size_t i = 0; i < waypoints.size(); ++i) {
            double dx = waypoints[i].x - current_x;
            double dy = waypoints[i].y - current_y;
            double dz = waypoints[i].z - current_z;
            waypoints[i].distance = std::sqrt(dx * dx + dy * dy + dz * dz);
        }
    }

    void sortByDistance() {
        std::sort(waypoints.begin(), waypoints.end(), compareByDistance);
    }

    void saveSortedWaypoints(const std::string& filename) {
        std::ofstream out(filename.c_str());
        if (!out.is_open()) return;

        out << std::fixed << std::setprecision(1);
        for (size_t i = 0; i < waypoints.size(); ++i) {
            out << waypoints[i].id << ","
                << waypoints[i].x << ","
                << waypoints[i].y << ","
                << waypoints[i].z << ","
                << waypoints[i].name << ","
                << waypoints[i].distance << "\n";
        }
    }

    void print() const {
        std::cout << std::fixed << std::setprecision(1);
        std::cout << "Sorted waypoints by distance:\n";
        for (size_t i = 0; i < waypoints.size(); ++i) {
            const Waypoint& w = waypoints[i];
            std::cout << "ID: " << w.id
                      << ", " << w.name
                      << ", Pos: (" << w.x << ", " << w.y << ", " << w.z << ")"
                      << ", Dist: " << w.distance << "\n";
        }
    }
};

int main() {
    WaypointSorter ws;

    if (!ws.loadWaypoints("waypoints.txt")) {
        std::cout << "Error: cannot open waypoints.txt\n";
        return 1;
    }

    double curX = 800.0, curY = 600.0, curZ = 1800.0;

    ws.calculateDistances(curX, curY, curZ);
    ws.sortByDistance();
    ws.saveSortedWaypoints("waypoints_sorted.txt");
    ws.print();

    std::cout << "Saved to waypoints_sorted.txt\n";
    return 0;
}
