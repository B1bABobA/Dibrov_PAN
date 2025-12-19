#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>
#include <fstream>
#include <cstdlib>

using namespace std;

const double MASS0 = 98000.0;
const double S_WING = 201.0;
const int ENGINE_COUNT = 3;
const double THRUST_PERCENT = 1;
const double TU154_NOMINAL_THRUST = 3 * 105000.0;

const double H_START = 300.0;
const double H_FINISH = 5000.0;
const double V_START_KMH = 300.0;
const double V_FINISH_KMH = 800.0;

const double DM = 0.25;
const double DH = 1000;
const double G = 9.81;
const double DEG_TO_RAD = 57.3;

const double MAX_VERTICAL_SPEED = 8.0;
const double MAX_CLIMB_ANGLE = 15.0;
const double MIN_CLIMB_SPEED_KMH = 350.0;

enum OptimizationCriterion {
    MIN_TIME = 1,
    MIN_FUEL = 2
};

enum ManeuverType {
    RAZGON = 1,
    PODIEM = 2,
    RAZGON_PODIEM = 3
};

struct AtmosPoint {
    double H, rho, a, T;
};

static const AtmosPoint ATMOS_TABLE[] = {
    {0.0,     1.22500, 340.294, 288.150},
    {500.0,   1.16727, 338.370, 284.900},
    {1000.0,  1.11166, 336.435, 281.651},
    {2000.0,  1.00655, 332.532, 275.154},
    {3000.0,  0.90254, 328.584, 268.659},
    {4000.0,  0.81935, 324.589, 262.166},
    {5000.0,  0.73643, 320.545, 255.676},
    {6000.0,  0.66011, 316.452, 249.187},
    {7000.0,  0.59002, 312.306, 242.700},
    {8000.0,  0.52678, 308.105, 236.215},
    {9000.0,  0.46706, 303.848, 229.733},
    {10000.0, 0.41351, 299.532, 223.252},
    {11000.0, 0.36480, 295.154, 216.774}
};
const int ATMOS_N = 13;

double interpolate(double x, double x0, double x1, double y0, double y1) {
    if (fabs(x1 - x0) < 1e-9) return y0;
    return y0 + (x - x0) * (y1 - y0) / (x1 - x0);
}

void atmosphere(double H, double& rho, double& a_sound) {
    if (H <= ATMOS_TABLE[0].H) {
        rho = ATMOS_TABLE[0].rho;
        a_sound = ATMOS_TABLE[0].a;
        return;
    }
    if (H >= ATMOS_TABLE[ATMOS_N - 1].H) {
        rho = ATMOS_TABLE[ATMOS_N - 1].rho;
        a_sound = ATMOS_TABLE[ATMOS_N - 1].a;
        return;
    }

    for (int i = 0; i < ATMOS_N - 1; i++) {
        if (H >= ATMOS_TABLE[i].H && H < ATMOS_TABLE[i + 1].H) {
            rho = interpolate(H, ATMOS_TABLE[i].H, ATMOS_TABLE[i + 1].H,
                ATMOS_TABLE[i].rho, ATMOS_TABLE[i + 1].rho);
            a_sound = interpolate(H, ATMOS_TABLE[i].H, ATMOS_TABLE[i + 1].H,
                ATMOS_TABLE[i].a, ATMOS_TABLE[i + 1].a);
            return;
        }
    }
}

bool is_in_flight_envelope(double H, double V_kmh) {
    if (V_kmh < 250.0 || V_kmh > 900.0) return false;
    if (H < 0.0 || H > 12000.0) return false;

    double rho, a_sound;
    atmosphere(H, rho, a_sound);
    double M = (V_kmh / 3.6) / a_sound;

    if (M > 0.85) return false;

    return true;
}

const double CY0 = -0.102;
const double CY1 = 0.085;

double Cx_alpha(double alpha_deg) {
    if (alpha_deg < 0.0) alpha_deg = 0.0;
    if (alpha_deg > 12.0) alpha_deg = 12.0;

    double Cx0 = 0.02;
    double k = 0.05;

    double Cy = CY0 + CY1 * alpha_deg;
    double Cx = Cx0 + k * Cy * Cy;

    if (Cx < 0.025) Cx = 0.025;

    return Cx;
}

double total_thrust(double H, double V_ms) {
    double thrust;
    return thrust = TU154_NOMINAL_THRUST * THRUST_PERCENT;
}

double specific_fuel_consumption(double H, double V_ms, double power_setting) {
    double rho, a_sound;
    atmosphere(H, rho, a_sound);
    double M = V_ms / a_sound;
    double H_km = H / 1000.0;

    double Cp_base = 0.70;

    double regime_factor;
    if (power_setting >= 1.0) {
        regime_factor = 1.0 + 0.50 * pow(power_setting - 1.0, 1.3);
    }
    else if (power_setting >= 0.88) {
        regime_factor = 0.93 - 0.03 * (power_setting - 0.88) / 0.12;
    }
    else if (power_setting >= 0.70) {
        regime_factor = 0.93 + 0.10 * pow((0.88 - power_setting) / 0.18, 1.2);
    }
    else {
        regime_factor = 1.15;
    }

    double altitude_factor = 1.0 - 0.08 * min(1.0, H_km / 11.0);
    double mach_factor = 1.0 + 0.15 * max(0.0, M - 0.5);

    double Cp = Cp_base * regime_factor * altitude_factor * mach_factor;
    return Cp / 9.81;
}

double calculate_alpha(double H, double V_ms, double mass) {
    double rho, a_sound;
    atmosphere(H, rho, a_sound);

    double P = total_thrust(H, V_ms);
    double q = 0.5 * rho * V_ms * V_ms;

    if (q < 100.0) return 8.0;

    double alpha_deg = (mass * G - P / DEG_TO_RAD - CY0 * q * S_WING) / (CY1 * q * S_WING);

    if (alpha_deg < 0.0) alpha_deg = 0.0;
    if (alpha_deg > 12.0) alpha_deg = 12.0;

    return alpha_deg;
}

struct SegmentData {
    double time;
    double fuel;
    bool valid;
};

SegmentData calculate_razgon(double H, double V1_ms, double V2_ms, double mass, double power_setting) {
    SegmentData result;
    result.valid = false;
    result.time = 1e9;
    result.fuel = 1e9;

    if (!is_in_flight_envelope(H, V1_ms * 3.6) || !is_in_flight_envelope(H, V2_ms * 3.6)) {
        return result;
    }

    double V_avg = 0.5 * (V1_ms + V2_ms);
    double alpha_deg = calculate_alpha(H, V_avg, mass);
    double alpha_rad = alpha_deg / DEG_TO_RAD;

    double P_max = total_thrust(H, V_avg);
    double P_used = P_max * power_setting;

    double rho, a_sound;
    atmosphere(H, rho, a_sound);
    double q = 0.5 * rho * V_avg * V_avg;

    double Cx = Cx_alpha(alpha_deg);
    double X = Cx * q * S_WING;

    double dV_dt = (P_used * cos(alpha_rad) - X) / mass;

    if (dV_dt <= 0.01) return result;

    double dt = (V2_ms - V1_ms) / dV_dt;

    if (dt > 1000.0 || dt <= 0) return result;

    double c_p = specific_fuel_consumption(H, V_avg, power_setting);
    double fuel = c_p * P_used * dt / 3600.0;

    result.time = dt;
    result.fuel = fuel;
    result.valid = true;

    return result;
}

SegmentData calculate_podiem(double H1, double H2, double V_ms, double mass, double power_setting, double max_vy_factor) {
    SegmentData result;
    result.valid = false;
    result.time = 1e9;
    result.fuel = 1e9;

    if (V_ms * 3.6 < MIN_CLIMB_SPEED_KMH) return result;
    if (!is_in_flight_envelope(H1, V_ms * 3.6) || !is_in_flight_envelope(H2, V_ms * 3.6)) {
        return result;
    }

    double H_avg = 0.5 * (H1 + H2);
    double alpha_deg = calculate_alpha(H_avg, V_ms, mass);

    double P_max = total_thrust(H_avg, V_ms);
    double P_used = P_max * power_setting;

    double rho, a_sound;
    atmosphere(H_avg, rho, a_sound);
    double q = 0.5 * rho * V_ms * V_ms;

    double Cx = Cx_alpha(alpha_deg);
    double X = Cx * q * S_WING;

    double P_excess = P_used - X;

    if (P_excess <= 0) return result;

    double theta_max_rad = MAX_CLIMB_ANGLE / DEG_TO_RAD;
    double sin_theta = min(P_excess / (mass * G), sin(theta_max_rad));

    if (sin_theta <= 0.005) return result;

    double Vy = V_ms * sin_theta;

    double max_vy_limit = MAX_VERTICAL_SPEED * max_vy_factor;
    if (Vy > max_vy_limit) {
        Vy = max_vy_limit;
    }

    double dt = (H2 - H1) / Vy;

    if (dt <= 0 || dt > 2000.0) return result;

    double c_p = specific_fuel_consumption(H_avg, V_ms, power_setting);
    double fuel = c_p * P_used * dt / 3600.0;

    result.time = dt;
    result.fuel = fuel;
    result.valid = true;

    return result;
}

SegmentData calculate_razgon_podiem(double H1, double H2, double V1_ms, double V2_ms, double mass, double power_setting, double max_vy_factor) {
    SegmentData result;
    result.valid = false;
    result.time = 1e9;
    result.fuel = 1e9;

    double V_avg = 0.5 * (V1_ms + V2_ms);
    double H_avg = 0.5 * (H1 + H2);

    if (V_avg * 3.6 < (MIN_CLIMB_SPEED_KMH * 0.95)) return result;
    if (!is_in_flight_envelope(H1, V1_ms * 3.6) || !is_in_flight_envelope(H2, V2_ms * 3.6)) {
        return result;
    }

    SegmentData razgon = calculate_razgon(H_avg, V1_ms, V2_ms, mass, power_setting);
    SegmentData podiem = calculate_podiem(H1, H2, V_avg, mass, power_setting, max_vy_factor);

    if (!razgon.valid || !podiem.valid) return result;

    double dH = H2 - H1;
    double dV_kmh = (V2_ms - V1_ms) * 3.6;

    double time_for_climb = dH / 5.0;
    double time_for_accel = fabs(dV_kmh) / 15.0;

    double dt = max(time_for_climb, time_for_accel);

    if (dt <= 0 || dt > 3000.0) return result;

    double Vy = dH / dt;
    double max_vy_limit = MAX_VERTICAL_SPEED * max_vy_factor * 1.2;

    if (Vy < 0.3 || Vy > max_vy_limit) return result;

    double dV_dt = (V2_ms - V1_ms) / dt;
    if (fabs(dV_dt) > 5.0) return result;

    double P_max = total_thrust(H_avg, V_avg);
    double P_used = P_max * power_setting;

    double c_p = specific_fuel_consumption(H_avg, V_avg, power_setting);
    double fuel = c_p * P_used * dt / 3600.0;

    result.time = dt;
    result.fuel = fuel;
    result.valid = true;

    return result;
}

struct TrajectoryResult {
    vector<pair<double, double>> path;
    vector<ManeuverType> maneuvers;
    double total_time;
    double total_fuel;
    double avg_vy;
    int used_razgon;
    int used_podiem;
    int used_combined;
};

TrajectoryResult solve_trajectory(OptimizationCriterion criterion) {
    TrajectoryResult trajectory;

    cout << "\n========================================\n";
    if (criterion == MIN_TIME) {
        cout << "KRITERII: MINIMIZACIA VREMENI\n";
    }
    else {
        cout << "KRITERII: MINIMIZACIA TOPLIVA\n";
    }
    cout << "========================================\n\n";

    double dH = DH;
    double dV_kmh = 1152.0 / 4.0;
    double NH_steps = (H_FINISH - H_START) / dH+1;
    double NV_steps = (V_FINISH_KMH - V_START_KMH) / dV_kmh+1;

    // Количество точек в сетке
    int NH_points = static_cast<int>(NH_steps) + 1;
    int NV_points = static_cast<int>(NV_steps) + 1;

    vector<double> H_grid(NH_points);
    vector<double> V_grid_kmh(NV_points);
    vector<double> V_grid_ms(NV_points);

    for (int i = 0; i < NH_points; i++) {
        H_grid[i] = H_START + i * dH;
    }
    for (int i = 0; i < NV_points; i++) {
        V_grid_kmh[i] = V_START_KMH + i * dV_kmh;
        V_grid_ms[i] = V_grid_kmh[i] / 3.6;
    }

    vector<double> power_settings;
    double max_vy_factor;

    if (criterion == MIN_TIME) {
        power_settings.push_back(1.08);
        power_settings.push_back(1.05);
        power_settings.push_back(1.00);
        max_vy_factor = 1.0;
    }
    else {
        power_settings.push_back(0.88);
        power_settings.push_back(0.82);
        power_settings.push_back(0.75);
        max_vy_factor = 0.65;
    }

    vector<vector<double> > cost_table(NH_points, vector<double>(NV_points, 1e9));
    vector<vector<double> > time_table(NH_points, vector<double>(NV_points, 0.0));
    vector<vector<double> > fuel_table(NH_points, vector<double>(NV_points, 0.0));
    vector<vector<int> > prev_i(NH_points, vector<int>(NV_points, -1));
    vector<vector<int> > prev_j(NH_points, vector<int>(NV_points, -1));
    vector<vector<ManeuverType> > maneuver_type(NH_points, vector<ManeuverType>(NV_points, RAZGON));

    cost_table[0][0] = 0.0;

    for (int i = 0; i < NH_points; i++) {
        for (int j = 0; j < NV_points; j++) {
            if (cost_table[i][j] >= 1e9) continue;

            double H1 = H_grid[i];
            double V1_ms = V_grid_ms[j];

            for (size_t ps = 0; ps < power_settings.size(); ps++) {
                double power_setting = power_settings[ps];

                if (j < NV_points) {
                    double V2_ms = V_grid_ms[j + 1];
                    SegmentData seg = calculate_razgon(H1, V1_ms, V2_ms, MASS0, power_setting);

                    if (seg.valid) {
                        double cost_increment = (criterion == MIN_TIME) ? seg.time : seg.fuel;
                        double new_cost = cost_table[i][j] + cost_increment;

                        if (new_cost < cost_table[i][j + 1]) {
                            cost_table[i][j + 1] = new_cost;
                            time_table[i][j + 1] = time_table[i][j] + seg.time;
                            fuel_table[i][j + 1] = fuel_table[i][j] + seg.fuel;
                            prev_i[i][j + 1] = i;
                            prev_j[i][j + 1] = j;
                            maneuver_type[i][j + 1] = RAZGON;
                        }
                    }
                }

                if (i < NH_points) {
                    double H2 = H_grid[i + 1];
                    SegmentData seg = calculate_podiem(H1, H2, V1_ms, MASS0, power_setting, max_vy_factor);

                    if (seg.valid) {
                        double cost_increment = (criterion == MIN_TIME) ? seg.time : seg.fuel;
                        double new_cost = cost_table[i][j] + cost_increment;

                        if (new_cost < cost_table[i + 1][j]) {
                            cost_table[i + 1][j] = new_cost;
                            time_table[i + 1][j] = time_table[i][j] + seg.time;
                            fuel_table[i + 1][j] = fuel_table[i][j] + seg.fuel;
                            prev_i[i + 1][j] = i;
                            prev_j[i + 1][j] = j;
                            maneuver_type[i + 1][j] = PODIEM;
                        }
                    }
                }

                if (i < NH_points - 1 && j < NV_points) {
                    double H2 = H_grid[i + 1];
                    double V2_ms = V_grid_ms[j + 1];
                    SegmentData seg = calculate_razgon_podiem(H1, H2, V1_ms, V2_ms, MASS0, power_setting, max_vy_factor);

                    if (seg.valid) {
                        double cost_increment = (criterion == MIN_TIME) ? seg.time : seg.fuel;
                        double new_cost = cost_table[i][j] + cost_increment;

                        if (new_cost < cost_table[i + 1][j + 1]) {
                            cost_table[i + 1][j + 1] = new_cost;
                            time_table[i + 1][j + 1] = time_table[i][j] + seg.time;
                            fuel_table[i + 1][j + 1] = fuel_table[i][j] + seg.fuel;
                            prev_i[i + 1][j + 1] = i;
                            prev_j[i + 1][j + 1] = j;
                            maneuver_type[i + 1][j + 1] = RAZGON_PODIEM;
                        }
                    }
                }
            }
        }
    }

    // Вывод матрицы времени
    cout << "Matrica vremeni (s):\n";
    cout << "     V->";
    for (int j = 0; j < NV_points; j++) {
        cout << setw(7) << (int)V_grid_kmh[j];
    }
    cout << "\nH\n";
    for (int i = 0; i < NH_points; i++) {
        cout << setw(5) << (int)H_grid[i];
        for (int j = 0; j < NV_points; j++) {
            if (time_table[i][j] < 1e8) {
                cout << setw(7) << (int)time_table[i][j];
            }
            else {
                cout << setw(7) << "---";
            }
        }
        cout << "\n";
    }

    // Вывод матрицы топлива
    cout << "\nMatrica raskhoda topliva (kg):\n";
    cout << "     V->";
    for (int j = 0; j < NV_points; j++) {
        cout << setw(7) << (int)V_grid_kmh[j];
    }
    cout << "\nH\n";
    for (int i = 0; i < NH_points; i++) {
        cout << setw(5) << (int)H_grid[i];
        for (int j = 0; j < NV_points; j++) {
            if (fuel_table[i][j] < 1e8) {
                cout << setw(7) << (int)fuel_table[i][j];
            }
            else {
                cout << setw(7) << "---";
            }
        }
        cout << "\n";
    }
    cout << "\n";

    // ИСПРАВЛЕНО: правильные индексы конечной точки
    int end_i = NH_points - 1;  // Последний индекс по высоте
    int end_j = NV_points - 1;  // Последний индекс по скорости

    if (cost_table[end_i][end_j] >= 1e9) {
        cout << "OSHIBKA: Ne naiden put!\n";
        return trajectory;
    }

    vector<pair<double, double> > path;
    vector<ManeuverType> path_maneuvers;
    int ci = end_i, cj = end_j;

    while (ci >= 0 && cj >= 0) {
        path.push_back(make_pair(H_grid[ci], V_grid_kmh[cj]));
        path_maneuvers.push_back(maneuver_type[ci][cj]);
        int pi = prev_i[ci][cj];
        int pj = prev_j[ci][cj];
        if (pi == -1) break;
        ci = pi;
        cj = pj;
    }

    reverse(path.begin(), path.end());
    reverse(path_maneuvers.begin(), path_maneuvers.end());

    int used_razgon = 0, used_podiem = 0, used_razgon_podiem = 0;
    for (size_t k = 1; k < path_maneuvers.size(); k++) {
        if (path_maneuvers[k] == RAZGON) used_razgon++;
        else if (path_maneuvers[k] == PODIEM) used_podiem++;
        else if (path_maneuvers[k] == RAZGON_PODIEM) used_razgon_podiem++;
    }

    cout << "Optimalnaya traektoriya:\n";
    cout << "-------------------------------------------------------------\n";
    cout << "Tochka\tH (m)\t\tV (km/h)\tManevr\n";
    cout << "-------------------------------------------------------------\n";

    for (size_t k = 0; k < path.size(); k++) {
        cout << k + 1 << "\t" << fixed << setprecision(1)
            << setw(8) << path[k].first
            << "\t" << setw(8) << path[k].second;

        if (k > 0) {
            if (path_maneuvers[k] == RAZGON) cout << "\t[Razgon]";
            else if (path_maneuvers[k] == PODIEM) cout << "\t[Podiem]";
            else if (path_maneuvers[k] == RAZGON_PODIEM) cout << "\t[Razgon+Podiem]";
        }
        else {
            cout << "\t[Start]";
        }
        cout << "\n";
    }

    cout << "\n=============================================\n";
    cout << "Ispolzovano v traektorii:\n";
    cout << "- Razgon: " << used_razgon << " raz\n";
    cout << "- Podiem: " << used_podiem << " raz\n";
    cout << "- Razgon+Podiem: " << used_razgon_podiem << " raz\n";
    cout << "---------------------------------------------\n";
    cout << fixed << setprecision(2);
    cout << "Vremya manevra:     " << time_table[end_i][end_j] << " s  ("  // ИСПРАВЛЕНО
        << time_table[end_i][end_j] / 60.0 << " min)\n";
    cout << "Raskhod topliva:    " << fuel_table[end_i][end_j] << " kg\n";  // ИСПРАВЛЕНО

    double delta_H = H_FINISH - H_START;
    double avg_climb_rate = delta_H / time_table[end_i][end_j];  // ИСПРАВЛЕНО

    cout << "Srednyaya Vy:       " << avg_climb_rate << " m/s  ("
        << avg_climb_rate * 60.0 << " m/min)\n";
    cout << "=============================================\n";

    trajectory.path = path;
    trajectory.maneuvers = path_maneuvers;
    trajectory.total_time = time_table[end_i][end_j];  // ИСПРАВЛЕНО
    trajectory.total_fuel = fuel_table[end_i][end_j];  // ИСПРАВЛЕНО
    trajectory.avg_vy = avg_climb_rate;
    trajectory.used_razgon = used_razgon;
    trajectory.used_podiem = used_podiem;
    trajectory.used_combined = used_razgon_podiem;

    return trajectory;
}

void create_html_plot(const TrajectoryResult& traj_time, const TrajectoryResult& traj_fuel) {
    ofstream html("trajectory_plot.html");

    html << "<!DOCTYPE html>\n";
    html << "<html>\n<head>\n";
    html << "<meta charset='UTF-8'>\n";
    html << "<title>TU-154 Trajectory Optimization</title>\n";
    html << "<script src='https://cdn.jsdelivr.net/npm/chart.js'></script>\n";
    html << "<style>\n";
    html << "body { font-family: Arial, sans-serif; margin: 20px; background: #f5f5f5; }\n";
    html << ".container { max-width: 1400px; margin: 0 auto; background: white; padding: 20px; border-radius: 10px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }\n";
    html << "h1 { color: #2c3e50; text-align: center; }\n";
    html << ".stats { display: flex; justify-content: space-around; margin: 20px 0; }\n";
    html << ".stat-box { background: #ecf0f1; padding: 15px; border-radius: 5px; text-align: center; flex: 1; margin: 0 10px; }\n";
    html << ".stat-box h3 { margin: 0; color: #34495e; font-size: 14px; }\n";
    html << ".stat-box p { margin: 10px 0 0 0; font-size: 24px; font-weight: bold; color: #2980b9; }\n";
    html << ".chart-container { position: relative; height: 600px; margin: 20px 0; }\n";
    html << "</style>\n";
    html << "</head>\n<body>\n";

    html << "<div class='container'>\n";
    html << "<h1>&#x1F6E9;&#xFE0F; Optimizacija traektorii TU-154 (Variant #1)</h1>\n";

    html << "<div class='stats'>\n";

    html << "<div class='stat-box' style='background: #e3f2fd;'>\n";
    html << "<h3>&#x26A1; Minimizacija vremeni</h3>\n";
    html << "<p>" << fixed << setprecision(1) << traj_time.total_time / 60.0 << " min</p>\n";
    html << "<p style='font-size: 16px; color: #666;'>" << traj_time.total_fuel << " kg topliva</p>\n";
    html << "</div>\n";

    html << "<div class='stat-box' style='background: #e8f5e9;'>\n";
    html << "<h3>&#x1F4B0; Minimizacija topliva</h3>\n";
    html << "<p>" << fixed << setprecision(1) << traj_fuel.total_time / 60.0 << " min</p>\n";
    html << "<p style='font-size: 16px; color: #666;'>" << traj_fuel.total_fuel << " kg topliva</p>\n";
    html << "</div>\n";

    html << "<div class='stat-box' style='background: #fff3e0;'>\n";
    html << "<h3>&#x1F4CA; Ekonomija topliva</h3>\n";
    html << "<p>" << fixed << setprecision(1) << (1.0 - traj_fuel.total_fuel / traj_time.total_fuel) * 100 << "%</p>\n";
    html << "<p style='font-size: 16px; color: #666;'>+" << (traj_fuel.total_time - traj_time.total_time) / 60.0 << " min</p>\n";
    html << "</div>\n";

    html << "</div>\n";

    html << "<div class='chart-container'>\n";
    html << "<canvas id='trajectoryChart'></canvas>\n";
    html << "</div>\n";

    html << "<script>\n";
    html << "const ctx = document.getElementById('trajectoryChart').getContext('2d');\n";

    html << "const dataMinTime = [\n";
    for (size_t i = 0; i < traj_time.path.size(); i++) {
        html << "{x: " << traj_time.path[i].second << ", y: " << traj_time.path[i].first << "}";
        if (i < traj_time.path.size() - 1) html << ",\n";
    }
    html << "\n];\n";

    html << "const dataMinFuel = [\n";
    for (size_t i = 0; i < traj_fuel.path.size(); i++) {
        html << "{x: " << traj_fuel.path[i].second << ", y: " << traj_fuel.path[i].first << "}";
        if (i < traj_fuel.path.size() - 1) html << ",\n";
    }
    html << "\n];\n";

    html << "const chart = new Chart(ctx, {\n";
    html << "  type: 'line',\n";
    html << "  data: {\n";
    html << "    datasets: [\n";
    html << "      {\n";
    html << "        label: 'Minimizacija vremeni (" << fixed << setprecision(1) << traj_time.total_time / 60.0 << " min, " << traj_time.total_fuel << " kg)',\n";
    html << "        data: dataMinTime,\n";
    html << "        borderColor: 'rgb(255, 99, 132)',\n";
    html << "        backgroundColor: 'rgba(255, 99, 132, 0.1)',\n";
    html << "        borderWidth: 3,\n";
    html << "        pointRadius: 5,\n";
    html << "        pointBackgroundColor: 'rgb(255, 99, 132)',\n";
    html << "        pointBorderColor: '#fff',\n";
    html << "        pointBorderWidth: 2,\n";
    html << "        tension: 0.1\n";
    html << "      },\n";
    html << "      {\n";
    html << "        label: 'Minimizacija topliva (" << fixed << setprecision(1) << traj_fuel.total_time / 60.0 << " min, " << traj_fuel.total_fuel << " kg)',\n";
    html << "        data: dataMinFuel,\n";
    html << "        borderColor: 'rgb(75, 192, 192)',\n";
    html << "        backgroundColor: 'rgba(75, 192, 192, 0.1)',\n";
    html << "        borderWidth: 3,\n";
    html << "        pointRadius: 5,\n";
    html << "        pointBackgroundColor: 'rgb(75, 192, 192)',\n";
    html << "        pointBorderColor: '#fff',\n";
    html << "        pointBorderWidth: 2,\n";
    html << "        tension: 0.1\n";
    html << "      }\n";
    html << "    ]\n";
    html << "  },\n";
    html << "  options: {\n";
    html << "    responsive: true,\n";
    html << "    maintainAspectRatio: false,\n";
    html << "    plugins: {\n";
    html << "      title: {\n";
    html << "        display: true,\n";
    html << "        text: 'Traektorija nabora vysoty H(V)',\n";
    html << "        font: { size: 18 }\n";
    html << "      },\n";
    html << "      legend: {\n";
    html << "        display: true,\n";
    html << "        position: 'top',\n";
    html << "        labels: { font: { size: 14 } }\n";
    html << "      },\n";
    html << "      tooltip: {\n";
    html << "        callbacks: {\n";
    html << "          label: function(context) {\n";
    html << "            return context.dataset.label + ': H=' + context.parsed.y.toFixed(0) + ' m, V=' + context.parsed.x.toFixed(0) + ' km/h';\n";
    html << "          }\n";
    html << "        }\n";
    html << "      }\n";
    html << "    },\n";
    html << "    scales: {\n";
    html << "      x: {\n";
    html << "        type: 'linear',\n";
    html << "        title: { display: true, text: 'Skorost V (km/h)', font: { size: 16 } },\n";
    html << "        min: 300,\n";
    html << "        max: 900,\n";
    html << "        grid: { color: 'rgba(0, 0, 0, 0.1)' }\n";
    html << "      },\n";
    html << "      y: {\n";
    html << "        type: 'linear',\n";
    html << "        title: { display: true, text: 'Vysota H (m)', font: { size: 16 } },\n";
    html << "        min: 0,\n";
    html << "        max: 8000,\n";
    html << "        grid: { color: 'rgba(0, 0, 0, 0.1)' }\n";
    html << "      }\n";
    html << "    }\n";
    html << "  }\n";
    html << "});\n";
    html << "</script>\n";

    html << "</div>\n";
    html << "</body>\n</html>\n";

    html.close();

    cout << "\n==============================================\n";
    cout << "Grafik sokhranyen v 'trajectory_plot.html'\n";
    cout << "Otkryvayu v brauzere...\n";
    cout << "==============================================\n\n";

	#ifdef _WIN32
    	system("start trajectory_plot.html");
	#elif __APPLE__
    	system("open trajectory_plot.html");
	#else
    	system("xdg-open trajectory_plot.html");
	#endif
}

int main() {
    cout << fixed << setprecision(2);

    cout << "\n=================================================\n";
    cout << "   OPTIMIZACIA TRAEKTORII TU-154 (Variant #1)\n";
    cout << "=================================================\n";
    cout << "Samolyet: TU-154 (massa " << MASS0 / 1000.0 << " t)\n";
    cout << "Start: H = " << H_START << " m, V = " << V_START_KMH << " km/h\n";
    cout << "Finish: H = " << H_FINISH << " m, V = " << V_FINISH_KMH << " km/h\n";
    cout << "=================================================\n\n";

    TrajectoryResult traj_time = solve_trajectory(MIN_TIME);
    TrajectoryResult traj_fuel = solve_trajectory(MIN_FUEL);
    create_html_plot(traj_time, traj_fuel);

    return 0;
}
