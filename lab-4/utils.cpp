#include "utils.h"
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>

double f(double t, double x, double y, double z, double a, double b) {
    return a * x + y + 10 * y * z;
}

double g(double t, double x, double y, double z, double a, double b) {
    return -x - 0.4 * y + 5*x*z;
}

double h(double t, double x, double y, double z, double a, double b) {
    return b*z - 5 * x * y;
}

bool ab_conditions(double a, double b) {
    return a + b < 0.4;
}

bool is_corr(double x, double y, double z) {
    return !(isnan(x) || isnan(y) || isnan(z) || isinf(x) || isinf(y) || isinf(z));
}

std::vector<double> get_stat_points(double a, double b)
{
    std::vector<double> stat_points = { 0,0,0 };
    double sq = 225. - 80 * a;
    if (sq < 0)
        return stat_points;
    if (sq < 1e-8) {
        double z = 0.05;
        if (b > 0)
            return stat_points;
        if (b > -1e-8) {
            if (is_corr(0, 0, z)) {
                stat_points.push_back(0);
                stat_points.push_back(0);
                stat_points.push_back(z);
            }
            return stat_points;
        }
        sq = - 0.4 * b / 3.;
        double x_1 = 0.2 * sqrt(sq);
        double x_2 = -x_1;
        double y_1 = -a * x_1 / (1 + 10 * z);
        double y_2 = -a * x_2 / (1 + 10 * z);

        if (is_corr(x_1, y_1, z)) {
            stat_points.push_back(x_1);
            stat_points.push_back(y_1);
            stat_points.push_back(z);
        }
        if (is_corr(x_2, y_2, z)) {
            stat_points.push_back(x_2);
            stat_points.push_back(y_2);
            stat_points.push_back(z);
        }

        return stat_points;
    }
    double z_1 = 0.05 + sqrt(sq)/100;
    double z_2 = 0.05 - sqrt(sq) / 100;

    // -----------
    sq = 2 * b * z_1 / (5 * z_1 - 1);
    if (sq > 0) {
        double x_1 = 0.2 * sqrt(sq);
        double x_2 = -x_1;
        double y_1 = -a * x_1 / (1 + 10 * z_1);
        double y_2 = -a * x_2 / (1 + 10 * z_1);

        if (is_corr(x_1, y_1, z_1)) {
            stat_points.push_back(x_1);
            stat_points.push_back(y_1);
            stat_points.push_back(z_1);
        }

        if (is_corr(x_2, y_2, z_1)) {
            stat_points.push_back(x_2);
            stat_points.push_back(y_2);
            stat_points.push_back(z_1);
        }
    }
    else if (sq > -1e-8) {
        if (is_corr(0, 0, z_1)) {
            stat_points.push_back(0);
            stat_points.push_back(0);
            stat_points.push_back(z_1);
        }
    }
    // -----------
    sq = 2 * b * z_2 / (5 * z_2 - 1);
    if (sq > 0) {
        double x_1 = 0.2 * sqrt(sq);
        double x_2 = -x_1;
        double y_1 = -a * x_1 / (1 + 10 * z_2);
        double y_2 = -a * x_2 / (1 + 10 * z_2);
        if (is_corr(x_1, y_1, z_2)) {
            stat_points.push_back(x_1);
            stat_points.push_back(y_1);
            stat_points.push_back(z_2);
        }
        if (is_corr(x_2, y_2, z_2)) {
            stat_points.push_back(x_2);
            stat_points.push_back(y_2);
            stat_points.push_back(z_2);
        }
    }
    else if (sq > -1e-8) {
        if (is_corr(0, 0, z_1)) {
            stat_points.push_back(0);
            stat_points.push_back(0);
            stat_points.push_back(z_2);
        }
    }
    return stat_points;
}

std::vector<double> new_point(double f(double, double, double, double, double, double),
    double g(double, double, double, double, double, double),
    double h(double, double, double, double, double, double),
    double t, double x, double y, double z,
    double step, double a, double b)
{
    double kx0, ky0, kz0, kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3;

    kx0 = f(t, x, y, z, a, b);
    ky0 = g(t, x, y, z, a, b);
    kz0 = h(t, x, y, z, a, b);

    kx1 = f(t + (step / 2), x + (step / 2) * kx0, y + (step / 2) * ky0, z + (step / 2) * kz0, a, b);
    ky1 = g(t + (step / 2), x + (step / 2) * kx0, y + (step / 2) * ky0, z + (step / 2) * kz0, a, b);
    kz1 = h(t + (step / 2), x + (step / 2) * kx0, y + (step / 2) * ky0, z + (step / 2) * kz0, a, b);

    kx2 = f(t + (step / 2), x + (step / 2) * kx1, y + (step / 2) * ky1, z + (step / 2) * kz1, a, b);
    ky2 = g(t + (step / 2), x + (step / 2) * kx1, y + (step / 2) * ky1, z + (step / 2) * kz1, a, b);
    kz2 = h(t + (step / 2), x + (step / 2) * kx1, y + (step / 2) * ky1, z + (step / 2) * kz1, a, b);

    kx3 = f(t + step, x + step * kx2, y + step * ky2, z + step * kz2, a, b);
    ky3 = g(t + step, x + step * kx2, y + step * ky2, z + step * kz2, a, b);
    kz3 = h(t + step, x + step * kx2, y + step * ky2, z + step * kz2, a, b);

    x = x + (step / 6) * (kx0 + 2 * kx1 + 2 * kx2 + kx3);
    y = y + (step / 6) * (ky0 + 2 * ky1 + 2 * ky2 + ky3);
    z = z + (step / 6) * (kz0 + 2 * kz1 + 2 * kz2 + kz3);

    return std::vector<double>{x, y, z};
}

void write_file(double a, double b,
    std::vector<double> t,
    std::vector<double> x,
    std::vector<double> y,
    std::vector<double> z,
    std::vector<double> stat_points,
    std::string filename)
{
    std::ofstream fout;
    fout.open(filename);

    if (t.size() != x.size() || t.size() != y.size() || t.size() != z.size()) {
        std::cout << "Dimension t is not equal with other dimensions" << std::endl;
        fout.close();
        throw "dimensions error";
    }

    fout << a << " " << b << "\n";
    fout << (int)(stat_points.size() / 3) << "\n";
    for (size_t i = 0; i < (int)(stat_points.size() / 3); i++) {
        fout << stat_points[3 * i] << " " << stat_points[3 * i + 1] << " " << stat_points[3 * i + 2] << "\n";
    }
    for (int i = 0; i < x.size(); i++) {
        fout << t[i] << " " << x[i] << " " << y[i] << " " << z[i] << "\n";
    }

    fout.close();
}