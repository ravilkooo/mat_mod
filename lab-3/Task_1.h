#include <vector>

long double FPU_der(long double r);
void verle(std::vector<long double>& q_cur,
           std::vector<long double>& v_cur,
           long long k,
           long long n,
           long double dt);
void simpl_verle(std::vector<long double>& q_cur,
           std::vector<long double>& v_cur,
           long long k,
           long long n,
           long double dt);