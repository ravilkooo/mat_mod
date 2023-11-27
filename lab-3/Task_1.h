#include <vector>

void do_experiment(long double alpha_in, long double beta_in,
    long double k, long double algo_type,
    std::vector<long long> times,
    std::string name,
    bool console_log);
long double FPU_der(long double r);
long double FPU(long double r);
void verle(std::vector<long double>& q_cur,
    std::vector<long double>& v_cur,
    long long k,
    long long n,
    long double dt,
    std::vector<long long> times,
    std::vector<std::vector<long double>>& checkpoints,
    std::vector<long double>& energy,
    bool console_log);
void simpl_verle(std::vector<long double>& q_cur,
    std::vector<long double>& v_cur,
    long long k,
    long long n,
    long double dt,
    std::vector<long long> times,
    std::vector<std::vector<long double>>& checkpoints,
    std::vector<long double>& energy,
    bool console_log);