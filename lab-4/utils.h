#pragma once
#include <vector>
#include <string>

double f(double t, double x, double y, double z, double a, double b);
double g(double t, double x, double y, double z, double a, double b);
double h(double t, double x, double y, double z, double a, double b);

bool ab_conditions(double a, double b);

std::vector<double> get_stat_points(double a, double b);

std::vector<double> new_point(double f(double, double, double, double, double, double), double g(double, double, double, double, double, double), double h(double, double, double, double, double, double), double t, double x, double y, double z, double step, double a, double b);

void write_file(double a, double b, std::vector<double> t, std::vector<double> x, std::vector<double> y, std::vector<double> z,
    std::vector<double> stat_points, std::string filename);
