#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <cmath>

struct Waypoint {
    int id;
    double x;
    double y;
    double z;
    double speed;
    std::string desc;
};

class WaypointManager {
public:
    WaypointManager(const std::string& filename = "waypoints.txt")
        : filename_(filename), current_index_(0), reach_tol_(1.0) {}

    void addWaypoint(int id, double x, double y, double z, double speed, const std::string& desc) {
        Waypoint w;
        w.id = id;
        w.x = x;
        w.y = y;
        w.z = z;
        w.speed = speed;
        w.desc = desc;
        waypoints_.push_back(w);
    }

    bool saveRoute() {
        std::ofstream out(filename_.c_str());
        if (!out.is_open()) return false;

        out << std::fixed << std::setprecision(1);
        for (size_t i = 0; i < waypoints_.size(); ++i) {
            out << waypoints_[i].id << ","
                << waypoints_[i].x << ","
                << waypoints_[i].y << ","
                << waypoints_[i].z << ","
                << waypoints_[i].speed << ","
                << waypoints_[i].desc << "\n";
        }
        return true;
    }

    bool loadRoute() {
        waypoints_.clear();
        current_index_ = 0;

        std::ifstream in(filename_.c_str());
        if (!in.is_open()) return false;

        std::string line;
        while (std::getline(in, line)) {
            if (line.empty()) continue;

            std::stringstream ss(line);
            std::string token;
            std::vector<std::string> parts;

            while (std::getline(ss, token, ',')) {
                while (!token.empty() && (token[0] == ' ' || token[0] == '\t')) token.erase(token.begin());
                while (!token.empty() && (token[token.size() - 1] == ' ' || token[token.size() - 1] == '\t')) token.pop_back();
                if (!token.empty() && token[0] == '\\') token = token.substr(1);
                parts.push_back(token);
            }

            if (parts.size() != 6) continue;

            Waypoint w;
            w.id = std::atoi(parts[0].c_str());
            w.x = std::atof(parts[1].c_str());
            w.y = std::atof(parts[2].c_str());
            w.z = std::atof(parts[3].c_str());
            w.speed = std::atof(parts[4].c_str());
            w.desc = parts[5];

            waypoints_.push_back(w);
        }

        return !waypoints_.empty();
    }

    double calculateTotalDistance() {
        if (waypoints_.size() < 2) return 0.0;

        double sum = 0.0;
        for (size_t i = 1; i < waypoints_.size(); ++i) {
            double dx = waypoints_[i].x - waypoints_[i - 1].x;
            double dy = waypoints_[i].y - waypoints_[i - 1].y;
            double dz = waypoints_[i].z - waypoints_[i - 1].z;
            sum += std::sqrt(dx * dx + dy * dy + dz * dz);
        }
        return sum;
    }

    Waypoint getNextWaypoint() {
        if (waypoints_.empty()) return Waypoint{0, 0, 0, 0, 0, ""};
        if (current_index_ + 1 >= waypoints_.size()) return waypoints_[current_index_];
        return waypoints_[current_index_ + 1];
    }

    bool checkWaypointReached(double x, double y, double z) {
        if (waypoints_.empty()) return false;

        Waypoint cur = waypoints_[current_index_];
        double dx = x - cur.x;
        double dy = y - cur.y;
        double dz = z - cur.z;
        double d = std::sqrt(dx * dx + dy * dy + dz * dz);

        if (d <= reach_tol_) {
            if (current_index_ + 1 < waypoints_.size()) current_index_++;
            return true;
        }
        return false;
    }

    bool hasWaypoints() const {
        return !waypoints_.empty();
    }

    size_t count() const {
        return waypoints_.size();
    }

    Waypoint currentWaypoint() const {
        if (waypoints_.empty()) return Waypoint{0, 0, 0, 0, 0, ""};
        return waypoints_[current_index_];
    }

private:
    std::string filename_;
    std::vector<Waypoint> waypoints_;
    size_t current_index_;
    double reach_tol_;
};

int main() {
    WaypointManager wm("waypoints.txt");

    if (!wm.loadRoute()) {
        std::cout << "Load error\n";
        return 1;
    }

    double total = wm.calculateTotalDistance();
    Waypoint cur = wm.currentWaypoint();
    Waypoint next = wm.getNextWaypoint();

    std::cout << std::fixed << std::setprecision(1);
    std::cout << "Route loaded: " << wm.count() << " waypoints\n";
    std::cout << "Total distance: " << total << " meters\n";
    std::cout << "Current waypoint: " << cur.desc << " (" << cur.x << ", " << cur.y << ", " << cur.z << ")\n";
    if (wm.count() >= 2 && (next.id != cur.id || next.desc != cur.desc)) {
        std::cout << "Next waypoint: " << next.desc << " (" << next.x << ", " << next.y << ", " << next.z << ")\n";
    } else {
        std::cout << "Next waypoint: none\n";
    }

    return 0;
}
