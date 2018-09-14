#include "BackgroundSolver.hpp"

BackgroundSolution solve_equations(Potential* pot, double N_star, double lim)
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
            N_temp += (pot->V(dx * n) / pot->dV(dx * n)) * dx;
        N_temp += 0.5 * (pot->V(phi_p) / pot->dV(phi_p)) * dx;
    }
    
    double dphi_p = - pot->dV(phi_p) / (3 * sqrt(pot->V(phi_p) / 3));
    std::vector<double> x0 = {phi_p, dphi_p, 0};
    
    //Find N_end
    double t = t0;
    std::vector<double> x = x0;
    dlsodar desolver(3, 1, 10000);
    desolver.integrate(t, 1e10, &x[0], equations, inflation_end, static_cast<void*>(ptrs));
    std::vector<double> x_end = x;
    double N_end = x_end[2];
    
    //Find Scalar Extrema
    std::vector<double> _N_extrema;
    t = t0;
    x = x0;
    params[0] = N_end;
    desolver = dlsodar(3, 3, 10000);
    while(x[2] < N_end)
    {
        desolver.integrate(t, 1e10, &x[0], equations, Extrema_Scalar, static_cast<void*>(ptrs));
        if(x[2] < N_end)
            _N_extrema.push_back(x[2]);
    }
    
    //Find Tensor Extrema
    std::vector<double> _N_extrema_tensor;
    t = t0;
    x = x0;
    params[0] = N_end;
    desolver = dlsodar(3, 3, 100000);
    while(x[2] < N_end)
    {
        desolver.integrate(t, 1e10, &x[0], equations, Extrema_Tensor, static_cast<void*>(ptrs));
        if(x[2] < N_end)
            _N_extrema_tensor.push_back(x[2]);
    }
    
    LinearInterpolator<double, double> _omega_2 = Solve_Variable(t0, x0, omega_2, _N_extrema, ptrs, lim);
    
    LinearInterpolator<double, double> _omega_2_tensor = Solve_Variable(t0, x0, omega_2_tensor, _N_extrema_tensor, ptrs, lim);
    LinearInterpolator<double, double> _dphi_H = Solve_Variable(t0, x0, dphi_H, _N_extrema, ptrs, lim);
    LinearInterpolator<double, double> _log_aH = Solve_Variable(t0, x0, log_aH, _N_extrema, ptrs, lim);
    
    double aH_star = exp(_log_aH(N_end - N_star));
    
    return BackgroundSolution(_omega_2, _omega_2_tensor, _log_aH, _dphi_H, _N_extrema, _N_extrema_tensor, aH_star, N_end);
}

BackgroundSolution solve_equations(Potential* pot, double N_star, double N_dagger, double lim)
{
    void* ptrs[2];
    ptrs[0] = static_cast<void*> (pot);
    double params[2];
    ptrs[1] = static_cast<void*> (params);
    
    double t0 = 1;
    double dphi_p = - sqrt(2.0/3);
    
    double phi_p = 0, N_temp_d = 0, N_end = 0, NN = 0;
    while(N_end < N_star + N_dagger)
    {
        phi_p += 0.5;
        double t = t0;
        std::vector<double> x = {phi_p, dphi_p, 0};
        dlsodar desolver(3, 1, 1000000);
        desolver.integrate(t, 1e10, &x[0], equations, inflation_end, static_cast<void*> (ptrs));
        N_end = x[2];
    }
    double phi_old = phi_p - 0.5;
    phi_p += 0.001;
    while(abs(N_temp_d - N_dagger) > 0.01)
    {
        double t = t0;
        std::vector<double> x = {phi_p, dphi_p, 0};
        dlsodar desolver(3, 1, 1000000);
        desolver.integrate(t, 1e10, &x[0], equations, inflation_end, static_cast<void*> (ptrs));
        N_end = x[2];
        
        t = t0;
        x = {phi_p, dphi_p, 0};
        desolver = dlsodar(3, 1, 100000);
        desolver.integrate(t, 1e10, &x[0], equations, inflation_begin, static_cast<void*> (ptrs));
        NN = x[2];
        
        params[0] = N_end - N_star;
        ptrs[1] = static_cast<void*> (params);
        desolver.integrate(t, 1e10, &x[0], equations, Find_N, static_cast<void*> (ptrs));
        N_temp_d = x[2] - NN;
        
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
    
    std::vector<double> x0 = {phi_p, dphi_p, 0};
    
    //Find Scalar Extrema
    std::vector<double> _N_extrema;
    auto t = t0;
    auto x = x0;
    params[0] = N_end;
    dlsodar desolver(3, 3, 10000);
    while(x[2] < N_end)
    {
        desolver.integrate(t, 1e10, &x[0], equations, Extrema_Scalar, static_cast<void*>(ptrs));
        if(x[2] < N_end)
            _N_extrema.push_back(x[2]);
    }
    
    //Find Tensor Extrema
    std::vector<double> _N_extrema_tensor;
    t = t0;
    x = x0;
    params[0] = N_end;
    desolver = dlsodar(3, 3, 100000);
    while(x[2] < N_end)
    {
        desolver.integrate(t, 1e10, &x[0], equations, Extrema_Tensor, static_cast<void*>(ptrs));
        if(x[2] < N_end)
            _N_extrema_tensor.push_back(x[2]);
    }
    
    LinearInterpolator<double, double> _omega_2 = Solve_Variable(t0, x0, omega_2, _N_extrema, ptrs, lim);    
    LinearInterpolator<double, double> _omega_2_tensor = Solve_Variable(t0, x0, omega_2_tensor, _N_extrema_tensor, ptrs, lim);
    LinearInterpolator<double, double> _dphi_H = Solve_Variable(t0, x0, dphi_H, _N_extrema, ptrs, lim);
    LinearInterpolator<double, double> _log_aH = Solve_Variable(t0, x0, log_aH, _N_extrema, ptrs, lim);
    
    double aH_star = exp(_log_aH(N_end - N_star));
    
    return BackgroundSolution(_omega_2, _omega_2_tensor, _log_aH, _dphi_H, _N_extrema,_N_extrema_tensor, aH_star, N_end);
}


LinearInterpolator<double, double> Solve_Variable(double t0, std::vector<double> x0, std::function<double(const double x[], Potential* pot)> Var, std::vector<double> N_extrema, void* ptrs[], double lim)
{
    auto pot = static_cast<Potential*> (ptrs[0]);
    double params[2];
    ptrs[1] = static_cast<void*> (params);
    
    LinearInterpolator<double, double> _Var;
    
    //Find N_end
    double t = t0;
    std::vector<double> x = x0;
    double N_i = x[2];
    dlsodar desolver(3, 1, 1e5);
    desolver.integrate(t, 1e10, &x[0], equations, inflation_end, static_cast<void*>(ptrs));
    double N_f = x[2];
    
    //Initialize Extrema in a single pass.
    _Var[N_i] = NAN; _Var[N_f] = NAN;
    for (auto N : N_extrema) if (N>N_i and N<N_f) _Var[N] = NAN;
    
    auto iter = _Var.begin();
    iter->second = Var(&x0[0], pot);
    ++iter;
    
    desolver = dlsodar(3, 1, 1e5);
    t = t0;
    x = x0;
    while(iter != _Var.end())
    {
        params[0] = iter->first;
        desolver.integrate(t, 1e10, &x[0], equations, Find_N, static_cast<void*>(ptrs));
        iter->second = Var(&x[0], pot);
        ++iter;
    }

    //Fill as necessary
    iter = _Var.begin();
    desolver = dlsodar(3, 1);
    t = t0;
    x = x0;
    while(iter != std::prev(_Var.end()))
    {
        N_i =  iter->first;
        N_f =  std::next(iter)->first;
        x0 = x;
        t0 = t;

        auto N_m = (N_i + N_f) / 2;
        params[0] = N_m;
        desolver.integrate(t, 1e10, &x[0], equations, Find_N, static_cast<void*>(ptrs));
        auto x_m1 = x;

        auto Approx = _Var(N_m);
        auto True = Var(&x[0], pot);
        
        if(abs(True - Approx) < lim)
            ++iter;
        else
        {
            _Var[N_m] = True;
            x = x0;
            t = t0;
            desolver.reset();
        }
    }
    
    return _Var;
    
}

double log_aH(const double x[], Potential* pot)
{
    double phi = x[0], dphi = x[1], n = x[2];
    return n + log(H(x,pot));
}

double dlog_aH(const double x[], Potential* pot)
{
    double phi = x[0], dphi = x[1], n = x[2];
    return 1 -dphi*dphi/2/H(x,pot)/H(x,pot);
}

double dphi_H(const double x[], Potential* pot)
{
    double phi = x[0], dphi = x[1];
    return dphi / H(x,pot);
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

double omega_2_tensor(const double x[], Potential* pot)
{
    double phi = x[0], dphi = x[1];
    return (3 * pow(dphi/H(x, pot), 2) + 6) * (pow(dphi/H(x, pot), 2) - 6) / 16.0 - 0.5 * pot->dV(phi) * dphi / pow(H(x, pot), 3);
}

double d_omega_2_tensor(const double x[], Potential* pot)
{
    double phi = x[0], dphi = x[1];
    return (3 * pow(dphi , 6)) / (8. * pow(H(x, pot), 6)) - (3 * pow(dphi, 4)) / pow(H(x, pot), 4) + (9 * pow(dphi, 2))/(2. * pow(H(x, pot) ,2)) - (3 * pow(dphi, 3) * pot->dV(phi)) / (2. * pow(H(x, pot), 5)) + (3 * dphi * pot->dV(phi)) / pow(H(x, pot), 3) + pow(pot->dV(phi), 2) / (2. * pow(H(x, pot), 4)) - (pow(dphi, 2) * pot->ddV(phi)) / (2. * pow(H(x, pot), 4));
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

void Extrema_Scalar(double g[], const double, const double x[], void* data)
{
    auto ptr = static_cast<void**>(data);
    auto pot = static_cast<Potential*> (ptr[0]);
    auto params = static_cast<double*> (ptr[1]);
    
    g[0] = d_omega_2(&x[0], pot);
    g[1] = params[0] - x[2];
    g[2] = dlog_aH(&x[0],pot);
}

void Extrema_Tensor(double g[], const double, const double x[], void* data)
{
    auto ptr = static_cast<void**>(data);
    auto pot = static_cast<Potential*> (ptr[0]);
    auto params = static_cast<double*> (ptr[1]);
    
    g[0] = d_omega_2_tensor(&x[0], pot);
    g[1] = params[0] - x[2];
    g[2] = dlog_aH(&x[0],pot);
}
