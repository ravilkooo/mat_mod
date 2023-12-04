#include <iostream>
#include <stdio.h>
#include <vector>
#include <iomanip>
#include <time.h>
#include <string>
#include "utils.h"

int run(double x0, double y0, double z0, double a, double b, std::string filename);
int run_mesh();

int main() {
    run_mesh();
}

int run_mesh()
{
    int cnt = 0;
    //double a = -5, b = 45;
    double x0 = 1, y0 = 1, z0 = 1;
    double a_start = -2.5; double a_fin = 2.5; double a_step = 0.5;
    double b_start = -2.5; double b_fin = 2.5; double b_step = 0.5;
    for (int i = 0; i < (int) (a_fin - a_start) / a_step + 1; i++) {
        for (int j = 0; j < (int)(b_fin - b_start) / b_step + 1; j++) {
            std::string filename = "results3/result_" + std::to_string(i) + "_" + std::to_string(j);
            std::cout << "run " << i << "," << j << std::endl;
            std::cout << "\t(" << a_start + i * a_step << ", \t" << b_start + j * b_step << ")" << std::endl;
            if (run(x0, y0, z0, a_start + i * a_step, b_start + j * b_step, filename) > 1)
                cnt++;
        }
    }
    std::cout << "##################\n" << cnt << "\n##################\n";
    return cnt;
}

int run(double x0, double y0, double z0, double a, double b, std::string filename)
{
    if (!ab_conditions(a,b))
        return 0;
    double step = 0.001;
    int n = 100000;
    double t0 = 0, t1 = t0 + n * step;

    std::vector<double> t(n), x(n), y(n), z(n);
    std::vector<double> point(3);
    double t_;

    //std::cout << "Enter a, b, x0, y0, z0 : ";
    //std::cin >> a >> b >> x0 >> y0 >> z0;
    //std::cin >> x0 >> y0 >> z0;
    x[0] = x0; y[0] = y0; z[0] = z0;

    auto stat_points = get_stat_points(a, b);
    std::cout << (int) (stat_points.size() / 3) << " stat points\n";
    if ((int)(stat_points.size() / 3) < 2) {
        return 0;
    }
    /*for (size_t i = 0; i < (int)(stat_points.size() / 3); i++) {
        std::cout << stat_points[3 * i] << "\t" << stat_points[3 * i + 1] << "\t" << stat_points[3 * i + 2] << "\n";
    }*/

    for (int i = 0; i < n - 1; i++) {
        t_ = i * step;
        point = new_point(f, g, h, t_, x[i], y[i], z[i], step, a, b);
        t[i + 1] = t_;
        x[i + 1] = point[0];
        y[i + 1] = point[1];
        z[i + 1] = point[2];
    }

    write_file(a, b, t, x, y, z, stat_points, filename);
  
    return (int)(stat_points.size() / 3);
}
