#include "BackgroundSolver.hpp"

double log_aH(const double x[], Potential* pot)
{
    double phi = x[0], dphi = x[1], n = x[2];
    return n + log(sqrt((0.5 * dphi * dphi + pot->V(phi)) / 3.0));
}

double dphi_H(const double x[], Potential* pot)
{
    double phi = x[0], dphi = x[1];
    return dphi / sqrt((0.5 * dphi * dphi + pot->V(phi)) / 3.0);
}

double H(const double x[], Potential* pot)
{
    double phi = x[0], dphi = x[1];
    return sqrt((0.5 * dphi * dphi + pot->V(phi)) / 3.0);
}

double omega_2(const double x[], Potential* pot)
{
    double phi = x[0], dphi = x[1];
    return (pot->ddV(phi) + 1.5 * pot->dV(phi) * dphi / H(x, pot)) / pow(H(x, pot), 2) - (pow(dphi/H(x, pot), 2) - 6) * (5 * pow(dphi/H(x, pot), 2) - 6) / 16.0;
}

double d_omega_2(const double x[], Potential* pot)
{
    double phi = x[0], dphi = x[1];
    return -3 * pow(pot->dV(phi) / H(x, pot) / H(x, pot), 2) - 9 * pot->dV(phi) * dphi / pow(H(x, pot), 3) + (7.0/2) * pot->dV(phi) * pow(dphi, 3) / pow(H(x, pot), 5) + (5.0/2) * pot->ddV(phi) * pow(dphi / H(x, pot) / H(x, pot), 2) - (27.0/2) * pow(dphi / H(x, pot), 2) + 6 * pow(dphi / H(x, pot), 4) - (5.0/8) * pow(dphi / H(x, pot), 6) + pot->dddV(phi) * dphi / pow(H(x, pot), 3);
}

void equations(double dx_dt[], const double, const double x[], void* data)
{
    auto ptr = static_cast<void**>(data);
    auto pot = static_cast<Potential*> (ptr[0]);
    dx_dt[0] = x[1];
    dx_dt[1] = - (3 * H(x, pot) * x[1] + pot->dV(x[0]));
    dx_dt[2] = H(x, pot);
}


void check(double g[], const double, const double x[], void*)
{
    g[0] = x[0];
}

void inflation_end(double g[], const double, const double x[], void* data)
{
    auto ptr = static_cast<void**>(data);
    auto pot = static_cast<Potential*> (ptr[0]);
    auto params = static_cast<double*> (ptr[1]);
    
    auto p = -(2.5 * x[1] * x[1] + x[1] * pot->dV(x[0]) / H(x, pot) - pow(0.5 * x[1] * x[1] / H(x, pot), 2));
    if(p > 0)
        g[0] = (-H(x, pot) + 0.5 * x[1] * x[1] / H(x, pot));
    else if(p < 0)
        g[0] = p;
    
    //////Fix THIS//////
    g[1] = params[0] - x[2];
}

void inflation_begin(double g[], const double, const double x[], void* data)
{
    auto ptr = static_cast<void**>(data);
    auto pot = static_cast<Potential*> (ptr[0]);
    
    g[0] = (-H(x, pot) + 0.5 * x[1] * x[1] / H(x, pot));
    
}

void Find_N(double g[], const double, const double x[], void* data)
{
    auto ptr = static_cast<void**>(data);
    auto params = static_cast<double*> (ptr[1]);
    
    g[0] = params[0] - x[2];
}

BackgroundSolution solve_equations(Potential* pot, double N_star)
{
    void* ptrs[2];
    ptrs[0] = static_cast<void*> (pot);
    double params[2];
    ptrs[1] = static_cast<void*> (params);
    
    double t0 = 1;
    double phi_p = 0, N_temp = 0;
    while(abs(N_temp - (N_star + 20)) > 0.1)
    {
        phi_p += 0.01;
        int steps = 1000;
        auto dx = (phi_p - 1e-5) / steps;
        
        N_temp = 0.5 * (pot->V(1e-5) / pot->dV(1e-5)) * dx;
        for(int n = 1; n < steps; n++)
        {
            N_temp += (pot->V(dx * n) / pot->dV(dx * n)) * dx;
        }
        N_temp += 0.5 * (pot->V(phi_p) / pot->dV(phi_p)) * dx;
    }
    
    double dphi_p = - pot->dV(phi_p) / (3 * sqrt(pot->V(phi_p) / 3));
    std::vector<double> x0 = {23, -0.8, 0};
    
    //Find N_end
    double t = t0;
    std::vector<double> x = x0;
    dlsodar desolver(3, 1, 10000);
    desolver.integrate(t, 1e10, &x[0], equations, inflation_end, static_cast<void*>(ptrs));
    std::vector<double> x_end = x;
    double N_end = x_end[2];
    
    std::vector<double> _N_extrema;
    //Run Solver
    double dt = 1e-2;
    t = t0;
    x = x0;
    Solutions sol;
    sol(x, t);
    desolver = dlsodar(3, 1, 100000);
    double old = d_omega_2(&x[0], pot);
    while(x[2] < N_end)
    {
        desolver.integrate(t, t + dt, &x[0], equations, inflation_end, ptrs);
        sol(x, t);
        dt *= 1.001;
        double New = d_omega_2(&x[0], pot);
        if(old < 0 and New > 0)
        {
            _N_extrema.push_back(x[2]);
        }
        else if(old > 0 and New < 0)
        {
            _N_extrema.push_back(x[2]);
        }
        old = New;
    }

    LinearInterpolator<double, double> _omega_2 = Solve_Variable(t0, x0, omega_2, _N_extrema, ptrs, 1e-4);
    LinearInterpolator<double, double> _dphi_H = Solve_Variable(t0, x0, dphi_H, _N_extrema, ptrs, 1e-4);
    LinearInterpolator<double, double> _log_aH = Solve_Variable(t0, x0, log_aH, _N_extrema, ptrs, 1e-4);
    
    double aH_star = exp(_log_aH(N_end - N_star));
    
    return BackgroundSolution(_omega_2, _log_aH, _dphi_H, _N_extrema, aH_star, N_end);
}


LinearInterpolator<double, double> Solve_Variable(double t0, std::vector<double> x0, std::function<double(const double x[], Potential* pot)> Var, std::vector<double> N_extrema, void* ptrs[], double lim)
{
    auto pot = static_cast<Potential*> (ptrs[0]);
    double params[2];
    ptrs[1] = static_cast<void*> (params);
    
    std::vector<std::pair<double, double>> N_pair;
    LinearInterpolator<double, double> _Var;
    
    //Find N_end
    double t = t0;
    std::vector<double> x = x0;
    dlsodar desolver(3, 1, 10000);
    desolver.integrate(t, 1e10, &x[0], equations, inflation_end, static_cast<void*>(ptrs));
    std::vector<double> x_end = x;
    params[0] = x_end[2];
    
    double N_i = x0[2], N_f = x_end[2];
    
    //Initialize Extrema
    if(N_extrema.size() != 0)
    {
        auto p0 = static_cast<size_t>(std::lower_bound(N_extrema.begin(), N_extrema.end(), N_i) - N_extrema.begin());
        auto p1 = static_cast<size_t>(std::lower_bound(N_extrema.begin(), N_extrema.end(), N_f) - N_extrema.begin());
        if(p1 <= N_extrema.size())
        {
            N_pair.push_back(std::make_pair(N_i, N_extrema[p0]));
            for(size_t n = p0+1; n < p1; n++)
                N_pair.push_back(std::make_pair(N_extrema[n-1], N_extrema[n]));
            
            N_pair.push_back(std::make_pair(N_extrema[p1-1], N_f));
        }
    }
    
    if(N_pair.size() == 0)
    {
        N_pair.push_back(std::make_pair(N_i, N_f));
    }
    
    
    _Var.insert(N_pair[0].first, Var(&x0[0], pot));
    
    desolver = dlsodar(3, 2, 100000);
    t = t0;
    x = x0;
    for(size_t n = 1; n < N_pair.size(); n++)
    {
        params[0] = N_pair[n].first;
        desolver.integrate(t, 1e10, &x[0], equations, Find_N, ptrs);
        
        _Var.insert(N_pair[n].first, Var(&x[0], pot));
    }
    
    params[0] = N_pair.back().second;
    desolver.integrate(t, 1e10, &x[0], equations, Find_N, ptrs);
    
    _Var.insert(N_pair.back().second, Var(&x[0], pot));
    
    int count = 0;
    while(N_pair.size() != 0)
    {
        desolver = dlsodar(3, 2, 100000);
        for(size_t n = 0; n < N_pair.size(); n++)
        {
            N_i = N_pair[n].first;
            N_f = N_pair[n].second;
            
            t = t0;
            x = x0;
            auto N_m1 = (2 * N_i + N_f) / 3.0;
            params[0] = N_m1;
            desolver.integrate(t, 1e10, &x[0], equations, Find_N, ptrs);
            std::vector<double> x_m1 = x;
            
            auto N_m2 = (N_i + 2 * N_f) / 3.0;
            params[0] = N_m2;
            desolver.integrate(t, 1e10, &x[0], equations, Find_N, ptrs);
            std::vector<double> x_m2 = x;
            
            count += 2;
            
            auto temp_approx1 = _Var(N_m1);
            auto temp_approx2 = _Var(N_m2);
            auto temp_true1 = Var(&x_m1[0], pot);
            auto temp_true2 = Var(&x_m2[0], pot);
            
            N_pair.erase(N_pair.begin() + static_cast<int>(n));
            
            if(abs(temp_true1 - temp_approx1) > lim and abs(temp_true2 - temp_approx2) > lim)
            {
                N_pair.insert(N_pair.begin() + static_cast<int>(n), std::make_pair(N_i, N_m1));
                N_pair.insert(N_pair.begin() + static_cast<int>(n) + 1, std::make_pair(N_m1, N_m2));
                N_pair.insert(N_pair.begin() + static_cast<int>(n) + 2, std::make_pair(N_m2, N_f));
                n += 2;
                _Var.insert(N_m1, (temp_true1));
                _Var.insert(N_m2, (temp_true2));
            }
            if(abs(temp_true1 - temp_approx1) > lim and abs(temp_true2 - temp_approx2) < lim)
            {
                N_pair.insert(N_pair.begin() + static_cast<int>(n), std::make_pair(N_i, N_m1));
                N_pair.insert(N_pair.begin() + static_cast<int>(n) + 1, std::make_pair(N_m1, N_m2));
                n += 1;
                _Var.insert(N_m1, (temp_true1));
            }
            if(abs(temp_true1 - temp_approx1) < lim and abs(temp_true2 - temp_approx2) > lim)
            {
                N_pair.insert(N_pair.begin() + static_cast<int>(n), std::make_pair(N_m1, N_m2));
                N_pair.insert(N_pair.begin() + static_cast<int>(n) + 1, std::make_pair(N_m2, N_f));
                n += 1;
                _Var.insert(N_m2, (temp_true2));
            }
        }
    }
    
    return _Var;
    
}
