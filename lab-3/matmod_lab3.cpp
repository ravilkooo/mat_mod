#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>

#include "Task_1.h"

#define VARIANT 11
#define SIMPLEX 1
#define COMMON 0

long double a = 1;
long double m = 1;

long double beta, alpha;

void task3();
void task4(long long t_start, long long t_steps, long long t_step_len, int alpha_steps);
void task4_borders(long long t_start, long long t_end, long long t_step_len, int alpha_steps);
void anim(long long t_start, long long t_steps);

int main()
{
    //task3();
    //task4((long long)1, (long long) (1e6 / 5e2), 5e2, 10);
    task4_borders((long long)1e6, (long long) 2e6, 1e2, 10);
    //anim(1, 10000);
}

void task3()
{
    std::vector<long long> times =
        // {0, 1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000};
    { 0, 1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000 };
    do_experiment(0.8, 0.3, 100000000, SIMPLEX, times, "1e8_SIMPLEX", true);
    do_experiment(0.8, 0.3, 100000000, COMMON, times, "1e8_COMMON", true);
}

void task4(long long t_start, long long t_steps, long long t_step_len, int alpha_steps)
{
    std::vector<long long> times(t_steps+1, 0);
    for (int i = 0; i < times.size()-1; i++) {
        times[i + 1] = t_start + i * t_step_len;
    }
    for (int i = 0; i < alpha_steps-1; i++) {
        long double _alpha = (i + 1.) / alpha_steps;
        std::string name = "./task4/" + std::to_string(i+1);
        do_experiment(_alpha, 0., t_start+t_steps * t_step_len, SIMPLEX, times, name, false);
    }
}

void task4_borders(long long t_start, long long t_end, long long t_step_len, int alpha_steps)
{
    int t_steps = (int)((t_end - t_start) / t_step_len);
    std::vector<long long> times(t_steps + 1, 0);
    for (int i = 0; i < times.size() - 1; i++) {
        times[i + 1] = t_start + i * t_step_len;
    }
    for (int i = 0; i < alpha_steps - 1; i++) {
        long double _alpha = (i + 1.) / alpha_steps;
        std::string name = "./task4/" + std::to_string(i + 1);
        do_experiment(_alpha, 0., t_end+1, SIMPLEX, times, name, false);
    }
}

void anim(long long t_start, long long t_steps)
{
    std::vector<long long> times(t_steps + 1, 0);
    for (int i = 0; i < times.size() - 1; i++) {
        times[i + 1] = t_start + i*10;
    }
    do_experiment(0.5, 0., t_start + t_steps*10, SIMPLEX, times, "anima", false);
}

void do_experiment(long double alpha_in, long double beta_in,
    long double k, long double algo_type,
    std::vector<long long> times,
    std::string name,
    bool console_log)
{
    // k - timesteps
    alpha = alpha_in;
    beta = beta_in;

    std::cout << "START" << std::endl;
    std::cout << alpha << " " << beta << " " << k << " ";
    std::cout << ((algo_type == COMMON) ? "Common" : "Simplex") << std::endl;

    int n = 100;
    std::vector<std::vector<long double>> checkpoints(2 * times.size() + 2, std::vector<long double>(2 * n + 2));
    std::vector<long double> energy(times.size(), 0);

    long double dt = 0.01;
    long double t_end = k * dt;
    long double q = 0.01 * VARIANT;
    long double m = 1; // mass
    long double a = 1;

    // init values
    std::vector<long double> q_init(2 * n + 2, 0.0);
    std::vector<long double> v_init(2 * n + 2, 0.0);

    q_init[n] -= q;
    q_init[n + 1] += q;

    std::vector<long double> q_cur = std::vector<long double>(q_init);
    std::vector<long double> v_cur = std::vector<long double>(v_init);

    // calculation
    std::cout << "START CALCULATION" << std::endl;
    if (times[0] == 0)
    {
        checkpoints[0] = q_cur;
        checkpoints[1] = v_cur;

        for (int j = 1; j < 2 * n + 1; j++) {
            energy[0] += FPU(q_cur[j] - q_cur[j - 1] );
            energy[0] += 0.5 * v_cur[j] * v_cur[j];
        }
        if (console_log)
        {
            std::cout << "checkpoint: " << std::setw(10) << times[0] << " || ";
            std::cout << "E = " << std::setprecision(17) << energy[0] << std::endl;
        }
    }

    if (algo_type == COMMON)
        verle(q_cur, v_cur, k, n, dt, times, checkpoints, energy, console_log);
    else
        simpl_verle(q_cur, v_cur, k, n, dt, times, checkpoints, energy, console_log);
    std::cout << "FINISH CALCULATION" << std::endl;

    // saving
    //std::string dir = "results/";
    std::string code = name;
    if (code == "")
        code = (algo_type == COMMON ? "COM_" : "SMP_") + std::to_string((long long)(k + 100 * alpha + beta));

    std::string out_filename = code + "_pars" + ".out";
    //std::ofstream f_out(dir + out_filename);
    std::ofstream f_out(out_filename);
    f_out << k << " " << alpha << " " << beta << std::endl;
    f_out.close();

    out_filename = code + "_res" + ".out";
    //f_out.open(dir + out_filename);
    f_out.open(out_filename);
    f_out << "x";
    for (int i = 0; i < times.size(); i++) {
        f_out << ",q_" << i << ",v_" << i;
    }
    f_out << std::endl;
    for (int j = 1; j < 2 * n + 1; j++) {
        f_out << a * j;
        for (int i = 0; i < times.size(); i++) {
            f_out << "," << std::setw(21) << std::setprecision(17) << checkpoints[2 * i][j] << "," << std::setprecision(17) << checkpoints[2 * i + 1][j];
            // "," << q_cur[j] << "," << v_cur[j] << std::endl;
        }
        f_out << std::endl;
    }
    f_out.close();

    out_filename = code + "_energy" + ".out";
    //f_out.open(dir + out_filename);
    f_out.open(out_filename);
    for (int i = 0; i < times.size(); i++) {
        f_out << std::setprecision(17) << energy[i] << std::endl;
    }
    f_out << std::endl;
    f_out.close();

    std::cout << "FINISH" << std::endl;
    // free(x_sc_arr);
    // free(v_sc_arr);
}

int old_main()
{
    std::cout << "START" << std::endl;
    // std::cout << std::setprecision (17) << xi << std::endl;

    long long k = 1000000; // time
    int n = 100;
    int algo_type = SIMPLEX;

    std::cout << "Enter alpha and beta: ";
    std::cin >> alpha >> beta;
    std::cout << "Enter k (time_steps): ";
    std::cin >> k;
    std::cout << "Choose algorythm type: ";
    std::cout << "Simplex(" << SIMPLEX << ") or ";
    std::cout << "Common(" << COMMON << "): ";
    std::cin >> algo_type;

    long double dt = 0.01;
    long double t_end = k * dt;
    long double q = 0.01 * VARIANT;
    long double m = 1; // mass
    long double a = 1;

    // init values
    std::vector<long double> q_init(2 * n + 2, 0.0);
    std::vector<long double> v_init(2 * n + 2, 0.0);

    q_init[n] -= q;
    q_init[n + 1] += q;

    std::vector<long double> q_cur = std::vector<long double>(q_init);
    std::vector<long double> v_cur = std::vector<long double>(v_init);

    // calculation
    std::cout << "START CALCULATION" << std::endl;
    std::vector<long long> times;
    std::vector<std::vector<long double>> checkpoints;
    std::vector<long double> energy;
    if (algo_type == COMMON)
        verle(q_cur, v_cur, k, n, dt, times, checkpoints, energy, true);
    else
        simpl_verle(q_cur, v_cur, k, n, dt, times, checkpoints, energy, true);
    std::cout << "FINISH CALCULATION" << std::endl;

    // saving
    std::string dir = "results/";
    std::string code = (algo_type == COMMON ? "COM_" : "SMP_") + std::to_string((long long)(k + 100 * alpha + beta));

    std::string out_filename = code + "_pars" + ".out";
    std::ofstream f_out(dir + out_filename);
    f_out << k << " " << alpha << " " << beta << std::endl;
    f_out.close();

    out_filename = code + "_res" + ".out";
    f_out.open(dir + out_filename);
    f_out << "x,q,v" << std::endl;
    for (int i = 1; i < 2 * n + 1; i++) {
        f_out << a * i << "," << q_cur[i] << "," << v_cur[i] << std::endl;
    }
    f_out.close();

    std::cout << "FINISH" << std::endl;
    // free(x_sc_arr);
    // free(v_sc_arr);
    return 0;
}

long double FPU(long double r)
{
    //long double r_abs = abs(r);
    //return r_abs * (1 + r_abs * (alpha + beta * r_abs));
    return r * r * (0.5 + r * (alpha / 3 + beta * r * 0.25));
}

long double FPU_der(long double r)
{
    //long double r_abs = abs(r);
    //return r_abs * (1 + r_abs * (alpha + beta * r_abs));
    return r * (1 + r * (alpha + beta * r));
}

void verle(std::vector<long double>& q_cur,
    std::vector<long double>& v_cur,
    long long k,
    long long n,
    long double dt,
    std::vector<long long> times,
    std::vector<std::vector<long double>>& checkpoints,
    std::vector<long double>& energy,
    bool console_log)
{
    int curr_checkp = times[0] == 0 ? 1 : 0;

    std::vector<long double> a_cur = std::vector<long double>(v_cur);
    for (int i = 1; i < 2 * n + 1; i++) {
        a_cur[i] = -(-FPU_der(q_cur[i + 1] - q_cur[i] ) + FPU_der(q_cur[i] - q_cur[i - 1] ));
    }
    a_cur[0] = a_cur[2 * n];
    a_cur[2 * n + 1] = a_cur[1];

    std::vector<long double> v_half = std::vector<long double>(v_cur);
    for (long long i = 0; i < k; i++) {

        for (int j = 0; j < 2 * n + 2; j++) {
            q_cur[j] = q_cur[j] + dt * (v_cur[j] + 0.5 * a_cur[j] * dt);
        }

        for (int j = 0; j < 2 * n + 2; j++) {
            v_cur[j] = v_cur[j] + 0.5 * a_cur[j] * dt;
        }

        for (int j = 1; j < 2 * n + 1; j++) {
            a_cur[j] = -(-FPU_der(q_cur[j + 1] - q_cur[j] ) + FPU_der(q_cur[j] - q_cur[j - 1] ));
        }

        a_cur[0] = a_cur[2 * n];
        a_cur[2 * n + 1] = a_cur[1];

        for (int j = 0; j < 2 * n + 2; j++) {
            v_cur[j] = v_cur[j] + 0.5 * a_cur[j] * dt;
        }

        if (curr_checkp < times.size() && (i + 1 == times[curr_checkp])) {
            checkpoints[2 * curr_checkp] = q_cur;
            checkpoints[2 * curr_checkp + 1] = v_cur;

            for (int j = 1; j < 2 * n + 1; j++) {
                energy[curr_checkp] += FPU(q_cur[j] - q_cur[j - 1] );
                energy[curr_checkp] += 0.5 * v_cur[j] * v_cur[j];
            }

            if (console_log)
            {
                std::cout << "checkpoint: " << std::setw(10) << times[curr_checkp] << " || ";
                std::cout << "E = " << std::setw(21) << std::setprecision(17) << energy[curr_checkp] << " || ";
                std::cout << "dE = " << std::setw(21) << std::setprecision(17) << abs(energy[curr_checkp] - energy[0]) << std::endl;
            }
            curr_checkp += 1;
        }
    }
}

void simpl_verle(std::vector<long double>& q_cur,
    std::vector<long double>& v_cur,
    long long k,
    long long n,
    long double dt,
    std::vector<long long> times,
    std::vector<std::vector<long double>>& checkpoints,
    std::vector<long double>& energy,
    bool console_log)
{
    int curr_checkp = times[0] == 0 ? 1 : 0;

    long double xi = 0.1931833275037836;

    std::vector<long double> a_cur = std::vector<long double>(v_cur);
    for (int i = 1; i < 2 * n + 1; i++) {
        a_cur[i] = -(-FPU_der(q_cur[i + 1] - q_cur[i] ) + FPU_der(q_cur[i] - q_cur[i - 1] ));
    }
    a_cur[0] = a_cur[2 * n];
    a_cur[2 * n + 1] = a_cur[1];

    std::vector<long double> q_1 = std::vector<long double>(v_cur);
    std::vector<long double> v_1 = std::vector<long double>(v_cur);
    std::vector<long double> a_1 = std::vector<long double>(v_cur);

    std::vector<long double> q_2 = std::vector<long double>(v_cur);
    std::vector<long double> a_2 = std::vector<long double>(v_cur);

    for (long long i = 0; i < k; i++) {
        // 1
        for (int j = 0; j < 2 * n + 2; j++) {
            q_cur[j] = q_cur[j] + xi * v_cur[j] * dt;
        }
        // 2
        for (int j = 1; j < 2 * n + 1; j++) {
            a_cur[j] = -(-FPU_der(q_cur[j + 1] - q_cur[j] ) + FPU_der(q_cur[j] - q_cur[j - 1] ));
        }
        a_cur[0] = a_cur[2 * n];
        a_cur[2 * n + 1] = a_cur[1];
        // 3
        for (int j = 0; j < 2 * n + 2; j++) {
            v_cur[j] = v_cur[j] + 0.5 * a_cur[j] * dt;
        }
        // 4
        for (int j = 0; j < 2 * n + 2; j++) {
            q_cur[j] = q_cur[j] + (1 - 2 * xi) * v_cur[j] * dt;
        }
        // 5
        for (int j = 1; j < 2 * n + 1; j++) {
            a_cur[j] = -(-FPU_der(q_cur[j + 1] - q_cur[j] ) + FPU_der(q_cur[j] - q_cur[j - 1] ));
        }
        a_cur[0] = a_cur[2 * n];
        a_cur[2 * n + 1] = a_cur[1];
        // 6
        for (int j = 0; j < 2 * n + 2; j++) {
            v_cur[j] = v_cur[j] + 0.5 * a_cur[j] * dt;
        }
        // 7
        for (int j = 0; j < 2 * n + 2; j++) {
            q_cur[j] = q_cur[j] + xi * v_cur[j] * dt;
        }
        if (curr_checkp < times.size() && (i+1 == times[curr_checkp])) {
            checkpoints[2 * curr_checkp] = q_cur;
            checkpoints[2 * curr_checkp + 1] = v_cur;

            for (int j = 1; j < 2 * n + 1; j++) {
                energy[curr_checkp] += FPU(q_cur[j] - q_cur[j - 1] );
                energy[curr_checkp] += 0.5 * v_cur[j] * v_cur[j];
            }
            if (console_log)
            {
                std::cout << "checkpoint: " << std::setw(10) << times[curr_checkp] << " || ";
                std::cout << "E = " << std::setw(21) << std::setprecision(17) << energy[curr_checkp] << " || ";
                std::cout << "dE = " << std::setw(21) << std::setprecision(17) << abs(energy[curr_checkp] - energy[0]) << std::endl;
            }
            curr_checkp += 1;
        }
    }
}