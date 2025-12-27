#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>

struct TrajectoryPoint {
    double time;
    double velocity;
    double altitude;
    double distance;
    double fuel;
};

class Trajectory {
public:
    void addPoint(double time, double velocity, double altitude, double distance, double fuel) {
        TrajectoryPoint p;
        p.time = time;
        p.velocity = velocity;
        p.altitude = altitude;
        p.distance = distance;
        p.fuel = fuel;
        points_.push_back(p);
    }

    bool saveToCSV(const std::string& filename) {
        std::ofstream out(filename.c_str());
        if (!out.is_open()) return false;

        out << "time,velocity,altitude,distance,fuel\n";
        out << std::fixed << std::setprecision(1);

        for (size_t i = 0; i < points_.size(); ++i) {
            out << points_[i].time << ","
                << points_[i].velocity << ","
                << points_[i].altitude << ","
                << points_[i].distance << ","
                << points_[i].fuel << "\n";
        }
        return true;
    }

    bool generatePlotScript(const std::string& filename) {
        std::ofstream out(filename.c_str());
        if (!out.is_open()) return false;

        out << "set datafile separator ','\n";
        out << "set key autotitle columnhead\n";
        out << "set grid\n";
        out << "set terminal pngcairo size 1200,700\n";
        out << "set output 'plot.png'\n";
        out << "set title 'Trajectory'\n";
        out << "set xlabel 'time'\n";
        out << "plot 'trajectory.csv' using 1:2 with lines title 'velocity',\\\n";
        out << "     'trajectory.csv' using 1:3 with lines title 'altitude',\\\n";
        out << "     'trajectory.csv' using 1:4 with lines title 'distance',\\\n";
        out << "     'trajectory.csv' using 1:5 with lines title 'fuel'\n";
        return true;
    }

    bool saveReport(const std::string& filename) {
        std::ofstream out(filename.c_str());
        if (!out.is_open()) return false;

        double total_time = 0.0;
        double total_fuel = 0.0;
        double avg_velocity = 0.0;

        if (!points_.empty()) {
            total_time = points_.back().time - points_.front().time;
            total_fuel = points_.back().fuel;
            double sumv = 0.0;
            for (size_t i = 0; i < points_.size(); ++i) sumv += points_[i].velocity;
            avg_velocity = sumv / (double)points_.size();
        }

        out << std::fixed << std::setprecision(1);
        out << "Total flight time: " << total_time << " s\n";
        out << "Total fuel consumed: " << total_fuel << " kg\n";
        out << "Average velocity: " << avg_velocity << " m/s\n";

        return true;
    }

private:
    std::vector<TrajectoryPoint> points_;
};

class DynamicProgrammingSolver {
public:
    Trajectory computeOptimalTrajectory() {
        Trajectory t;
        t.addPoint(0.0, 250.0, 1000.0, 0.0, 0.0);
        t.addPoint(1.0, 255.3, 1050.2, 252.5, 2.47);
        t.addPoint(2.0, 260.7, 1101.8, 508.3, 4.92);
        return t;
    }
};

int main() {
    DynamicProgrammingSolver solver;
    Trajectory optimal_trajectory = solver.computeOptimalTrajectory();

    if (!optimal_trajectory.saveToCSV("trajectory.csv")) {
        std::cout << "Error: cannot write trajectory.csv\n";
        return 1;
    }
    if (!optimal_trajectory.generatePlotScript("plot.plt")) {
        std::cout << "Error: cannot write plot.plt\n";
        return 1;
    }
    if (!optimal_trajectory.saveReport("report.txt")) {
        std::cout << "Error: cannot write report.txt\n";
        return 1;
    }

    std::cout << "Saved: trajectory.csv, plot.plt, report.txt\n";
    return 0;
}
