#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>

struct Rule {
    std::string field;
    double minv;
    double maxv;
};

class DataValidator {
public:
    DataValidator(const std::string& report_file = "validation_report.txt")
        : report_file_(report_file) {}

    void addValidationRule(const std::string& field, double min, double max) {
        Rule r;
        r.field = field;
        r.minv = min;
        r.maxv = max;
        rules_.push_back(r);
    }

    bool validateCoordinates(double x, double y, double z) {
        last_x_ = x;
        last_y_ = y;
        last_z_ = z;
        coord_ok_ = true;

        const Rule* rx = findRule("x");
        const Rule* ry = findRule("y");
        const Rule* rz = findRule("z");

        if (rx && (x < rx->minv || x > rx->maxv)) coord_ok_ = false;
        if (ry && (y < ry->minv || y > ry->maxv)) coord_ok_ = false;
        if (rz && (z < rz->minv || z > rz->maxv)) coord_ok_ = false;

        if (rz && (z < rz->minv || z > rz->maxv)) {
            coord_reason_ = "Altitude " + fmt1(z) + " exceeds maximum " + fmt1(rz->maxv);
        } else if (rx && (x < rx->minv || x > rx->maxv)) {
            coord_reason_ = "X " + fmt1(x) + " out of range [" + fmt1(rx->minv) + ", " + fmt1(rx->maxv) + "]";
        } else if (ry && (y < ry->minv || y > ry->maxv)) {
            coord_reason_ = "Y " + fmt1(y) + " out of range [" + fmt1(ry->minv) + ", " + fmt1(ry->maxv) + "]";
        } else {
            coord_reason_ = "";
        }

        return coord_ok_;
    }

    bool validateSpeed(double speed) {
        last_speed_ = speed;
        speed_ok_ = true;

        const Rule* r = findRule("speed");
        if (r && (speed < r->minv || speed > r->maxv)) speed_ok_ = false;

        if (!speed_ok_ && r) {
            speed_reason_ = "Speed " + fmt1(speed) + " exceeds maximum " + fmt1(r->maxv);
        } else {
            speed_reason_ = "";
        }

        return speed_ok_;
    }

    bool validateAcceleration(double acceleration) {
        last_accel_ = acceleration;
        accel_ok_ = true;

        const Rule* r = findRule("acceleration");
        if (r && (acceleration < r->minv || acceleration > r->maxv)) accel_ok_ = false;

        if (!accel_ok_ && r) {
            accel_reason_ = "Acceleration " + fmt1(acceleration) + " exceeds maximum " + fmt1(r->maxv);
        } else {
            accel_reason_ = "";
        }

        return accel_ok_;
    }

    void generateValidationReport() {
        std::ofstream out(report_file_.c_str());
        if (!out.is_open()) return;

        out << std::fixed << std::setprecision(1);

        out << "Validation report:\n";

        if (coord_ok_) {
            out << "Coordinates: OK\n";
        } else {
            const Rule* rz = findRule("z");
            if (rz && (last_z_ < rz->minv || last_z_ > rz->maxv)) {
                out << "Coordinates: ERROR - altitude " << last_z_ << " exceeds maximum " << rz->maxv << "\n";
            } else if (!coord_reason_.empty()) {
                out << "Coordinates: ERROR - " << coord_reason_ << "\n";
            } else {
                out << "Coordinates: ERROR\n";
            }
        }

        if (speed_ok_) {
            out << "Speed: OK\n";
        } else {
            const Rule* r = findRule("speed");
            if (r) out << "Speed: ERROR - speed " << last_speed_ << " exceeds maximum " << r->maxv << "\n";
            else out << "Speed: ERROR\n";
        }

        if (accel_ok_) {
            out << "Acceleration: OK\n";
        } else {
            const Rule* r = findRule("acceleration");
            if (r) out << "Acceleration: ERROR - acceleration " << last_accel_ << " exceeds maximum " << r->maxv << "\n";
            else out << "Acceleration: ERROR\n";
        }

        out << "Overall result: " << std::fixed << std::setprecision(0) << getValidationScore()
            << "% of data valid\n";
    }

    double getValidationScore() {
        int ok = 0;
        int total = 3;
        if (coord_ok_) ok++;
        if (speed_ok_) ok++;
        if (accel_ok_) ok++;
        return (double)ok * 100.0 / (double)total;
    }

private:
    const Rule* findRule(const std::string& field) const {
        for (size_t i = 0; i < rules_.size(); ++i) {
            if (rules_[i].field == field) return &rules_[i];
        }
        return 0;
    }

    static std::string fmt1(double v) {
        std::ostringstream os;
        os << std::fixed << std::setprecision(1) << v;
        return os.str();
    }

private:
    std::string report_file_;
    std::vector<Rule> rules_;

    double last_x_ = 0.0, last_y_ = 0.0, last_z_ = 0.0;
    double last_speed_ = 0.0;
    double last_accel_ = 0.0;

    bool coord_ok_ = false;
    bool speed_ok_ = false;
    bool accel_ok_ = false;

    std::string coord_reason_;
    std::string speed_reason_;
    std::string accel_reason_;
};

int main() {
    DataValidator v("validation_report.txt");

    v.addValidationRule("x", -100000.0, 100000.0);
    v.addValidationRule("y", -100000.0, 100000.0);
    v.addValidationRule("z", -1000.0, 5000.0);
    v.addValidationRule("speed", 0.0, 300.0);
    v.addValidationRule("acceleration", 0.0, 20.0);

    double x = 5000.5, y = 3000.2, z = 10000.0;
    double speed = 350.0;
    double accel = 25.0;

    v.validateCoordinates(x, y, z);
    v.validateSpeed(speed);
    v.validateAcceleration(accel);

    v.generateValidationReport();

    std::cout << "Report saved to validation_report.txt\n";
    return 0;
}
