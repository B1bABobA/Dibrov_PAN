#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <sstream>

struct TelemetryData {
    double time;
    double altitude;
    double speed;
    double heading;
    double fuel;
};

class TelemetryLogger {
public:
    TelemetryLogger() : base_name_("telemetry_"), index_(1), max_records_(1000), current_records_(0) {
        openCurrentFile();
    }

    ~TelemetryLogger() {
        if (out_.is_open()) out_.close();
    }

    bool logData(double time, double altitude, double speed, double heading, double fuel) {
        if (!out_.is_open()) openCurrentFile();
        rotateFileIfNeeded();
        if (!out_.is_open()) return false;

        TelemetryData d;
        d.time = time;
        d.altitude = altitude;
        d.speed = speed;
        d.heading = heading;
        d.fuel = fuel;

        out_.write(reinterpret_cast<const char*>(&d), sizeof(TelemetryData));
        if (!out_) return false;

        current_records_++;
        total_records_++;
        last_ = d;
        if (total_records_ == 1) first_ = d;

        if (total_records_ == 1) {
            min_alt_ = max_alt_ = d.altitude;
            min_speed_ = max_speed_ = d.speed;
            min_fuel_ = max_fuel_ = d.fuel;
        } else {
            if (d.altitude < min_alt_) min_alt_ = d.altitude;
            if (d.altitude > max_alt_) max_alt_ = d.altitude;
            if (d.speed < min_speed_) min_speed_ = d.speed;
            if (d.speed > max_speed_) max_speed_ = d.speed;
            if (d.fuel < min_fuel_) min_fuel_ = d.fuel;
            if (d.fuel > max_fuel_) max_fuel_ = d.fuel;
        }

        sum_alt_ += d.altitude;
        sum_speed_ += d.speed;
        sum_heading_ += d.heading;
        sum_fuel_ += d.fuel;

        return true;
    }

    void rotateFileIfNeeded() {
        if (current_records_ < max_records_) return;
        if (out_.is_open()) out_.close();
        index_++;
        current_records_ = 0;
        openCurrentFile();
    }

    std::vector<TelemetryData> readLogFile(const std::string& filename) {
        std::vector<TelemetryData> res;
        std::ifstream in(filename.c_str(), std::ios::binary);
        if (!in.is_open()) return res;

        TelemetryData d;
        while (in.read(reinterpret_cast<char*>(&d), sizeof(TelemetryData))) {
            res.push_back(d);
        }
        return res;
    }

    void printLogSummary() {
        std::cout << "Telemetry log summary\n";
        std::cout << "Total records: " << total_records_ << "\n";
        if (total_records_ == 0) return;

        double avg_alt = sum_alt_ / (double)total_records_;
        double avg_speed = sum_speed_ / (double)total_records_;
        double avg_heading = sum_heading_ / (double)total_records_;
        double avg_fuel = sum_fuel_ / (double)total_records_;

        std::cout << std::fixed << std::setprecision(1);

        std::cout << "Time range: " << first_.time << " .. " << last_.time << "\n";
        std::cout << "Altitude min/max/avg: " << min_alt_ << " / " << max_alt_ << " / " << avg_alt << "\n";
        std::cout << "Speed    min/max/avg: " << min_speed_ << " / " << max_speed_ << " / " << avg_speed << "\n";
        std::cout << "Heading  avg: " << avg_heading << "\n";
        std::cout << "Fuel     min/max/avg: " << min_fuel_ << " / " << max_fuel_ << " / " << avg_fuel << "\n";
        std::cout << "Current file: " << currentFilename() << "\n";
    }

    std::string currentFilename() const {
        return buildFilename(index_);
    }

private:
    void openCurrentFile() {
        std::string fname = buildFilename(index_);
        out_.open(fname.c_str(), std::ios::binary | std::ios::app);
    }

    std::string buildFilename(int idx) const {
        std::ostringstream os;
        os << base_name_ << std::setw(3) << std::setfill('0') << idx << ".bin";
        return os.str();
    }

private:
    std::string base_name_;
    int index_;
    std::ofstream out_;
    size_t max_records_;
    size_t current_records_;

    size_t total_records_ = 0;

    TelemetryData first_{};
    TelemetryData last_{};

    double min_alt_ = 0.0, max_alt_ = 0.0;
    double min_speed_ = 0.0, max_speed_ = 0.0;
    double min_fuel_ = 0.0, max_fuel_ = 0.0;

    double sum_alt_ = 0.0;
    double sum_speed_ = 0.0;
    double sum_heading_ = 0.0;
    double sum_fuel_ = 0.0;
};

int main() {
    TelemetryLogger logger;

    logger.logData(0.0, 100.0, 25.0, 45.0, 80.0);
    logger.logData(1.0, 105.0, 27.0, 46.0, 79.5);
    logger.logData(2.0, 110.0, 29.0, 47.0, 79.0);

    logger.printLogSummary();

    std::vector<TelemetryData> data = logger.readLogFile(logger.currentFilename());
    std::cout << "\nRead back from file: " << data.size() << " records\n";
    for (size_t i = 0; i < data.size(); ++i) {
        std::cout << std::fixed << std::setprecision(1)
                  << "Time: " << data[i].time
                  << ", Altitude: " << data[i].altitude
                  << ", Speed: " << data[i].speed
                  << ", Heading: " << data[i].heading
                  << ", Fuel: " << data[i].fuel << "\n";
    }

    return 0;
}
