#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <algorithm>

struct Target {
    int id;
    std::string name;
    double x;
    double y;
    double z;
    double priority;
    double distance;
};

static bool cmpByDistance(const Target& a, const Target& b) {
    return a.distance < b.distance;
}

class TargetManager {
public:
    TargetManager() : filename_("targets.txt") {}

    void addTarget(int id, const std::string& name, double x, double y, double z,
                   double priority, double distance) {
        Target t;
        t.id = id;
        t.name = name;
        t.x = x;
        t.y = y;
        t.z = z;
        t.priority = priority;
        t.distance = distance;
        targets_.push_back(t);
    }

    bool removeTarget(int target_id) {
        for (size_t i = 0; i < targets_.size(); ++i) {
            if (targets_[i].id == target_id) {
                targets_.erase(targets_.begin() + (long long)i);
                return true;
            }
        }
        return false;
    }

    void saveTargetsToFile() {
        std::ofstream out(filename_.c_str());
        if (!out.is_open()) return;

        out << std::fixed << std::setprecision(1);
        for (size_t i = 0; i < targets_.size(); ++i) {
            out << targets_[i].id << ","
                << targets_[i].name << ","
                << targets_[i].x << ","
                << targets_[i].y << ","
                << targets_[i].z << ","
                << targets_[i].priority << ","
                << targets_[i].distance << "\n";
        }
    }

    void loadTargetsFromFile() {
        targets_.clear();

        std::ifstream in(filename_.c_str());
        if (!in.is_open()) return;

        std::string line;
        while (std::getline(in, line)) {
            if (line.empty()) continue;

            std::vector<std::string> parts;
            std::stringstream ss(line);
            std::string token;

            while (std::getline(ss, token, ',')) {
                while (!token.empty() && (token[0] == ' ' || token[0] == '\t')) token.erase(token.begin());
                while (!token.empty() && (token[token.size() - 1] == ' ' || token[token.size() - 1] == '\t')) token.pop_back();
                parts.push_back(token);
            }

            if (parts.size() != 7) continue;

            Target t;
            t.id = std::atoi(parts[0].c_str());
            t.name = parts[1];
            t.x = std::atof(parts[2].c_str());
            t.y = std::atof(parts[3].c_str());
            t.z = std::atof(parts[4].c_str());
            t.priority = std::atof(parts[5].c_str());
            t.distance = std::atof(parts[6].c_str());

            targets_.push_back(t);
        }
    }

    std::vector<Target> getHighPriorityTargets(double min_priority) {
        std::vector<Target> res;
        for (size_t i = 0; i < targets_.size(); ++i) {
            if (targets_[i].priority >= min_priority) res.push_back(targets_[i]);
        }
        return res;
    }

    void sortByDistance() {
        std::sort(targets_.begin(), targets_.end(), cmpByDistance);
    }

    void printAll() const {
        std::cout << "Targets in system:\n";
        std::cout << std::fixed << std::setprecision(1);
        for (size_t i = 0; i < targets_.size(); ++i) {
            const Target& t = targets_[i];
            std::cout << "ID: " << t.id << ", " << t.name
                      << ", Position: (" << t.x << ", " << t.y << ", " << t.z << ")"
                      << ", Priority: " << t.priority
                      << ", Distance: " << t.distance << "\n";
        }
    }

private:
    std::string filename_;
    std::vector<Target> targets_;
};

int main() {
    TargetManager manager;
    manager.loadTargetsFromFile();

    manager.addTarget(4, "Name4", 200.0, 100.0, 20.0, 0.7, 1800.0);

    manager.printAll();
    manager.saveTargetsToFile();

    return 0;
}
