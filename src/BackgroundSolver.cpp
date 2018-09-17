#include "BackgroundSolver.hpp"
#include <limits>
#include "utils.hpp"

BackgroundSolution solve_equations(double lim, Potential* pot, double N_star, double N_dagger)
{
    BackgroundSolution Bsol;
    void* ptrs[2];
    ptrs[0] = static_cast<void*> (pot);
    double params[2];
    ptrs[1] = static_cast<void*> (params);
    auto N_tot = N_star + N_dagger, N_end=0.;

    auto f = [&pot, &ptrs, N_tot, &N_end](double phi_p) -> double
    {
        auto dphi_p = -sqrt(2./3);
        double t = 1.;
        std::vector<double> x = {phi_p, dphi_p, 0.};
        dlsodar desolver(3, 1, 1e5);
        desolver.integrate(t, std::numeric_limits<double>::max(), &x[0], equations_t, inflating, static_cast<void*> (ptrs));
        auto N_start = x[2];
        desolver.integrate(t, std::numeric_limits<double>::max(), &x[0], equations_t, inflating, static_cast<void*> (ptrs));
        N_end = x[2];
        return (N_end - N_start) - N_tot;
    };

    auto phi_p = find_root<double>(f,0.1,100.,lim);

    auto dphi_p = - sqrt(2 / (pot->V(phi_p) + 1./3));
    std::vector<double> x0 = {phi_p, dphi_p};

    //Find Scalar Extrema
    double n = 0;
    auto x = x0;
    dlsodar desolver(2, 2, 1e5);
    while(n < N_end)
    {
        desolver.integrate(n, std::numeric_limits<double>::max(), &x[0], equations_n, Extrema_Scalar, static_cast<void*>(ptrs));
        if(n < N_end) Bsol.N_extrema.push_back(n);
    }

    //Find Tensor Extrema
    n = 0;
    x = x0;
    desolver = dlsodar(2, 2, 1e5);
    while(n < N_end)
    {
        desolver.integrate(n, std::numeric_limits<double>::max(), &x[0], equations_n, Extrema_Tensor, static_cast<void*>(ptrs));
        if(n < N_end) Bsol.N_extrema_tensor.push_back(n);
    }

    Bsol.omega_2 = Solve_Variable(0, N_end, x0, omega_2, Bsol.N_extrema, ptrs, lim);
    std::cout << Bsol.omega_2.size() << std::endl;
    Bsol.omega_2_tensor = Solve_Variable(0, N_end, x0, omega_2_tensor, Bsol.N_extrema_tensor, ptrs, lim);
    std::cout << Bsol.omega_2_tensor.size() << std::endl;
    Bsol.dphi_H = Solve_Variable(0, N_end, x0, dphi_H, {}, ptrs, lim);
    std::cout << Bsol.dphi_H.size() << std::endl;
    Bsol.log_aH = Solve_Variable(0, N_end, x0, log_aH, {}, ptrs, lim);
    std::cout << Bsol.log_aH.size() << std::endl;

    Bsol.aH_star = exp(Bsol.log_aH(N_end - N_star));
    Bsol.N_end = N_end;

    return Bsol;
}

SemiLogInterpolator<double, double> Solve_Variable(double N_i, double N_f, std::vector<double> x0, std::function<double(const double n, const double x[], Potential* pot)> Var, std::vector<double> N_extrema, void* ptrs[], double lim)
{
    auto pot = static_cast<Potential*> (ptrs[0]);
    double params[2];
    ptrs[1] = static_cast<void*> (params);

    SemiLogInterpolator<double, double> _Var;

    //Initialize Extrema in a single pass.
    _Var[N_i] = {NAN,NAN}; _Var[N_f] = {NAN,NAN};
    for (auto N : N_extrema) if (N>N_i and N<N_f) _Var[N] = {NAN,NAN};

    auto desolver = dlsodar(2, 0, 1e5);
    auto n = N_i; 
    auto x = x0;
    for(auto &v : _Var)
    {
        desolver.integrate(n, v.first, &x[0], equations_n, nullptr, static_cast<void*>(ptrs));
        v.second.first = Var(n, &x[0], pot);
    }

    //Fill as necessary
    auto iter = _Var.begin();
    desolver = dlsodar(2, 0, 1e5);
    n = N_i; x = x0;
    while(iter != std::prev(_Var.end()))
    {
        N_i =  iter->first; N_f =  std::next(iter)->first;
        auto N_m = (N_i + N_f) / 2;

        x0 = x; N_i = n;
        desolver.integrate(n, N_m, &x[0], equations_n, nullptr, static_cast<void*>(ptrs));

        auto True = Var(n, &x[0], pot);

        if(_Var.insert(iter,n,True,lim)) ++iter;
        else
        {
            desolver.reset();
            x = x0; n = N_i;
        }
    }

    return _Var;
}

double dphi_H(const double n, const double x[], Potential* pot)
{
    double phi = x[0], dphi = x[1];
    return dphi;
}

double log_aH(const double n, const double x[], Potential* pot)
{
    double phi = x[0], dphi = x[1];
    return n + log(H(n, x, pot));
}

double dlog_aH(const double n, const double x[], Potential* pot)
{
    double phi = x[0], dphi = x[1];
    return 1 - 0.5 * dphi * dphi;
}

double H(const double n, const double x[], Potential* pot)
{
    double phi = x[0], dphi = x[1];
    return sqrt(pot->V(phi) / (3. - 0.5 * dphi * dphi));
}

double H_t(const double x[], Potential* pot)
{
    double phi = x[0], dphi = x[1];
    return sqrt((dphi*dphi/2 + pot->V(phi)) / 3);
}

double omega_2(const double n, const double x[], Potential* pot)
{
    double phi = x[0], dphi = x[1];
    return (pot->ddV(phi) + 1.5 * pot->dV(phi) * dphi) / pow(H(n, x, pot), 2) - (pow(dphi, 2) - 6) * (5 * pow(dphi, 2) - 6) / 16.0;
}

double omega_2_tensor(const double n, const double x[], Potential* pot)
{
    double phi = x[0], dphi = x[1];
    return (3 * pow(dphi, 2) + 6) * (pow(dphi, 2) - 6) / 16.0 - 0.5 * pot->dV(phi) * dphi / pow(H(n, x, pot), 2);
}

double d_omega_2_tensor(const double n, const double x[], Potential* pot)
{
    double phi = x[0], dphi = x[1];
    return (3 * pow(dphi , 6)) / 8. - (3 * pow(dphi, 4)) + (9 * pow(dphi, 2)) / 2. - (3 * pow(dphi, 3) * pot->dV(phi)) / (2. * pow(H(n, x, pot), 2)) + (3 * dphi * pot->dV(phi)) / pow(H(n, x, pot), 2) + pow(pot->dV(phi), 2) / (2. * pow(H(n, x, pot), 4)) - (pow(dphi, 2) * pot->ddV(phi)) / (2. * pow(H(n, x, pot), 2));
}

double d_omega_2(const double n, const double x[], Potential* pot)
{
    double phi = x[0], dphi = x[1];
    return -3 * pow(pot->dV(phi) / H(n, x, pot) / H(n, x, pot), 2) - 9 * pot->dV(phi) * dphi / pow(H(n, x, pot), 2) + (7.0/2) * pot->dV(phi) * pow(dphi, 3) / pow(H(n, x, pot), 2) + (5.0/2) * pot->ddV(phi) * pow(dphi / H(n, x, pot), 2) - (27.0/2) * pow(dphi, 2) + 6 * pow(dphi , 4) - (5.0/8) * pow(dphi , 6) + pot->dddV(phi) * dphi / pow(H(n, x, pot), 2);
}

void equations_n(double dx_dn[], const double n, const double x[], void* data)
{
    auto ptr = static_cast<void**>(data);
    auto pot = static_cast<Potential*> (ptr[0]);
    dx_dn[0] = x[1];
    dx_dn[1] = - ((3 - 0.5 * x[1] * x[1] - 0.5 * pot->dV(x[0]) * x[1] / pot->V(x[0])) * x[1] + 3 * pot->dV(x[0]) / pot->V(x[0]));
}

void equations_t(double dx_dt[], const double t, const double x[], void* data)
{
    auto ptr = static_cast<void**>(data);
    auto pot = static_cast<Potential*> (ptr[0]);
    dx_dt[0] = x[1];
    dx_dt[1] = -3*H_t(x,pot)*x[1] - pot->dV(x[0]);
    dx_dt[2] = H_t(x,pot);
}

void inflating(double g[], const double t, const double x[], void* data)
{
    auto ptr = static_cast<void**>(data);
    auto pot = static_cast<Potential*> (ptr[0]);

    g[0] = pot->V(x[0]) - x[1]*x[1];
}

void Extrema_Scalar(double g[], const double n, const double x[], void* data)
{
    auto ptr = static_cast<void**>(data);
    auto pot = static_cast<Potential*> (ptr[0]);

    g[0] = d_omega_2(n, x, pot);
    g[1] = dlog_aH(n, x, pot);
}

void Extrema_Tensor(double g[], const double n, const double x[], void* data)
{
    auto ptr = static_cast<void**>(data);
    auto pot = static_cast<Potential*> (ptr[0]);

    g[0] = d_omega_2_tensor(n, x, pot);
    g[1] = dlog_aH(n, x, pot);
}
