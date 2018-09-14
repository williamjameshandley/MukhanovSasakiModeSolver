#include "BackgroundSolver.hpp"

BackgroundSolution solve_equations(Potential* pot, double N_star, double lim)
{
    void* ptrs[2];
    ptrs[0] = static_cast<void*> (pot);
    double params[2];
    ptrs[1] = static_cast<void*> (params);
    
    double phi_p = 0, N_temp = 0;
    while(abs(N_temp - (N_star + 20)) > 0.1)
    {
        phi_p += 0.01;
        int steps = 1000;
        auto dx = (phi_p - 1e-5) / steps;
        
        N_temp = 0.5 * (pot->V(1e-5) / pot->dV(1e-5)) * dx;
        for(int n = 1; n < steps; n++)
            N_temp += (pot->V(dx * n) / pot->dV(dx * n)) * dx;
        N_temp += 0.5 * (pot->V(phi_p) / pot->dV(phi_p)) * dx;
    }
    
    double dphi_p = - pot->dV(phi_p) / sqrt(pow(pot->V(phi_p), 2) - pow(pot->dV(phi_p), 2));
    std::vector<double> x0 = {phi_p, dphi_p};

    //Find N_end
    double n = 0;
    std::vector<double> x = x0;
    dlsodar desolver(2, 1, 1e5);
    desolver.integrate(n, 1000, &x[0], equations, inflation_end, static_cast<void*>(ptrs));
    std::vector<double> x_end = x;
    double N_end = n;
    
    //Find Scalar Extrema
    std::vector<double> _N_extrema;
    n = 0;
    x = x0;
    desolver = dlsodar(2, 2, 1e5);
    while(n < N_end)
    {
        desolver.integrate(n, 1000, &x[0], equations, Extrema_Scalar, static_cast<void*>(ptrs));
        if(n < N_end)
        {
            _N_extrema.push_back(n);
        }
    }
    
    //Find Tensor Extrema
    std::vector<double> _N_extrema_tensor;
    n = 0;
    x = x0;
    desolver = dlsodar(2, 2, 1e5);
    while(n < N_end)
    {
        desolver.integrate(n, 1000, &x[0], equations, Extrema_Tensor, static_cast<void*>(ptrs));
        if(n < N_end)
            _N_extrema_tensor.push_back(n);
    }
    
    LinearInterpolator<double, double> _omega_2 = Solve_Variable(0, x0, omega_2, _N_extrema, ptrs, lim);
    LinearInterpolator<double, double> _omega_2_tensor = Solve_Variable(0, x0, omega_2_tensor, _N_extrema_tensor, ptrs, lim);
    LinearInterpolator<double, double> _dphi_H = Solve_Variable(0, x0, dphi_H, _N_extrema, ptrs, lim);
    LinearInterpolator<double, double> _log_aH = Solve_Variable(0, x0, log_aH, _N_extrema, ptrs, lim);
    
    double aH_star = exp(_log_aH(N_end - N_star));
    
    return BackgroundSolution(_omega_2, _omega_2_tensor, _log_aH, _dphi_H, _N_extrema,_N_extrema_tensor, aH_star, N_end);
}

BackgroundSolution solve_equations(Potential* pot, double N_star, double N_dagger, double lim)
{
    void* ptrs[2];
    ptrs[0] = static_cast<void*> (pot);
    double params[2];
    ptrs[1] = static_cast<void*> (params);
    
    double phi_p = 0, dphi_p, N_temp_d = 0, N_end = 0, NN = 0;
    
    while(N_end < N_star + N_dagger)
    {
        phi_p += 0.5;
        dphi_p = - sqrt(2 / (pot->V(phi_p) + 1./3));
        double n = 0;
        std::vector<double> x = {phi_p, dphi_p};
        dlsodar desolver(2, 1, 1e5);
        desolver.integrate(n, 1000, &x[0], equations, inflation_end, static_cast<void*> (ptrs));
        N_end = n;
    }
    double phi_old = phi_p - 0.5;
    phi_p += 0.001;
    while(abs(N_temp_d - N_dagger) > 0.01)
    {
        dphi_p = - sqrt(2 / (pot->V(phi_p) + 1./3));
        double n = 0;
        std::vector<double> x = {phi_p, dphi_p};
        dlsodar desolver(2, 1, 1e5);
        desolver.integrate(n, 1000, &x[0], equations, inflation_end, static_cast<void*> (ptrs));
        N_end = n;
        
        n = 0;
        x = {phi_p, dphi_p};
        desolver = dlsodar(2, 1, 1e5);
        desolver.integrate(n, 1000, &x[0], equations, inflation_begin, static_cast<void*> (ptrs));
        NN = n;
        
        desolver.integrate(n, N_end - N_star, &x[0], equations, inflation_end, static_cast<void*> (ptrs));
        N_temp_d = n - NN;
        
        if(N_temp_d < N_dagger)
        {
            phi_old = phi_p;
            phi_p += 0.001;
        }
        else if(N_temp_d > N_dagger)
        {
            phi_p -= (phi_p - phi_old) / phi_p;
        }
        
    }
    
    dphi_p = - sqrt(2 / (pot->V(phi_p) + 1./3));
    std::vector<double> x0 = {phi_p, dphi_p};
    
    //Find Scalar Extrema
    std::vector<double> _N_extrema;
    double n = 0;
    auto x = x0;
    dlsodar desolver(2, 2, 1e5);
    while(n < N_end)
    {
        desolver.integrate(n, 1000, &x[0], equations, Extrema_Scalar, static_cast<void*>(ptrs));
        if(n < N_end)
        {
            _N_extrema.push_back(n);
        }
    }
    
    //Find Tensor Extrema
    std::vector<double> _N_extrema_tensor;
    n = 0;
    x = x0;
    desolver = dlsodar(2, 2, 1e5);
    while(n < N_end)
    {
        desolver.integrate(n, 1000, &x[0], equations, Extrema_Tensor, static_cast<void*>(ptrs));
        if(n < N_end)
            _N_extrema_tensor.push_back(n);
    }
    
    LinearInterpolator<double, double> _omega_2 = Solve_Variable(0, x0, omega_2, _N_extrema, ptrs, lim);
    LinearInterpolator<double, double> _omega_2_tensor = Solve_Variable(0, x0, omega_2_tensor, _N_extrema_tensor, ptrs, lim);
    LinearInterpolator<double, double> _dphi_H = Solve_Variable(0, x0, dphi_H, _N_extrema, ptrs, lim);
    LinearInterpolator<double, double> _log_aH = Solve_Variable(0, x0, log_aH, _N_extrema, ptrs, lim);
    
    double aH_star = exp(_log_aH(N_end - N_star));
    
    return BackgroundSolution(_omega_2, _omega_2_tensor, _log_aH, _dphi_H, _N_extrema,_N_extrema_tensor, aH_star, N_end);
}

LinearInterpolator<double, double> Solve_Variable(double n0, std::vector<double> x0, std::function<double(const double n, const double x[], Potential* pot)> Var, std::vector<double> N_extrema, void* ptrs[], double lim)
{
    auto pot = static_cast<Potential*> (ptrs[0]);
    double params[2];
    ptrs[1] = static_cast<void*> (params);
    
    LinearInterpolator<double, double> _Var;
    
    //Find N_end
    double n = n0;
    double N_i = n0;
    std::vector<double> x = x0;
    dlsodar desolver(2, 1, 1e5);
    desolver.integrate(n, 1000, &x[0], equations, inflation_end, static_cast<void*>(ptrs));
    double N_f = n;
    
    //Initialize Extrema in a single pass.
    _Var[N_i] = NAN; _Var[N_f] = NAN;
    for (auto N : N_extrema) if (N>N_i and N<N_f) _Var[N] = NAN;
    
    auto iter = _Var.begin();
    iter->second = Var(N_i, &x0[0], pot);
    ++iter;
    
    desolver = dlsodar(2, 1, 1e5);
    n = n0;
    x = x0;
    while(iter != _Var.end())
    {
        desolver.integrate(n, iter->first, &x[0], equations, inflation_end, static_cast<void*>(ptrs));
        iter->second = Var(n, &x[0], pot);
        ++iter;
    }
    
    //Fill as necessary
    iter = _Var.begin();
    desolver = dlsodar(2, 1, 1e5);
    n = n0;
    x = x0;
    while(iter != std::prev(_Var.end()))
    {
        N_i =  iter->first;
        N_f =  std::next(iter)->first;
        x0 = x;
        n0 = n;
        
        auto N_m = (N_i + N_f) / 2;
        desolver.integrate(n, N_m, &x[0], equations, inflation_end, static_cast<void*>(ptrs));
        auto x_m1 = x;

        auto Approx = _Var(N_m);
        auto True = Var(n, &x[0], pot);
        
        if(abs(True - Approx) < lim)
            ++iter;
        else
        {
            _Var[N_m] = True;
            x = x0;
            n = n0;
            desolver.reset();
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

void equations(double dx_dt[], const double n, const double x[], void* data)
{
    auto ptr = static_cast<void**>(data);
    auto pot = static_cast<Potential*> (ptr[0]);
    dx_dt[0] = x[1];
    dx_dt[1] = - ((3 - 0.5 * x[1] * x[1] - 0.5 * pot->dV(x[0]) * x[1] / pot->V(x[0])) * x[1] + 3 * pot->dV(x[0]) / pot->V(x[0]));
}

void inflation_end(double g[], const double n, const double x[], void* data)
{
    auto ptr = static_cast<void**>(data);
    auto pot = static_cast<Potential*> (ptr[0]);
    
    auto p = -(3 - 0.5 * x[1] * x[1] - 0.5 * x[1] * pot->dV(x[0]) / pot->V(x[0])) * x[1] * x[1] - 3 * pot->dV(x[0]) * x[1] / pot->V(x[0]);
    if(p > 0)
        g[0] = -(1 - 0.5 * x[1] * x[1]);
    else if(p < 0)
        g[0] = p;
}

void inflation_begin(double g[], const double n, const double x[], void* data)
{
    g[0] = -(1 - 0.5 * x[1] * x[1]);
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
