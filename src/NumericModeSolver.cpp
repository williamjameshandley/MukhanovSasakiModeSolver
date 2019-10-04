#include "ModeSolver.hpp"

double _H(double phi, double dphi, Potential* pot)
{
    return std::sqrt((0.5 * dphi * dphi + pot->V(phi)) / 3.0);
}

double dz_z(double phi, double dphi, Potential* pot)
{
    return -2 * _H(phi, dphi, pot) - pot->dV(phi) / dphi + dphi * dphi / (2 * _H(phi, dphi, pot));
}

void __equations(double dx_dt[], const double, const double x[], void* data)
{
    auto ptr = static_cast<void**>(data);
    auto pot = static_cast<Potential*> (ptr[0]);
    auto params = static_cast<double*> (ptr[1]);
    double log_k = std::log(params[2]);
    
    dx_dt[0] = x[1];
    dx_dt[1] = - (3 * _H(x[0], x[1], pot) * x[1] + pot->dV(x[0]));
    dx_dt[2] = _H(x[0], x[1], pot);
    dx_dt[3] = std::exp(-x[2]);
}

void _equations(double dx_dt[], const double, const double x[], void* data)
{
    auto ptr = static_cast<void**>(data);
    auto pot = static_cast<Potential*> (ptr[0]);
    auto params = static_cast<double*> (ptr[1]);
    double log_k = std::log(params[2]);
    
    dx_dt[0] = x[1];
    dx_dt[1] = - (3 * _H(x[0], x[1], pot) * x[1] + pot->dV(x[0]));
    dx_dt[2] = _H(x[0], x[1], pot);
    dx_dt[3] = x[5];
    dx_dt[4] = x[6];
    dx_dt[5] = -(_H(x[0], x[1], pot) + 2.0 * dz_z(x[0], x[1], pot)) * x[5] - std::exp(2 * (log_k - x[2])) * x[3];
    dx_dt[6] = -(_H(x[0], x[1], pot) + 2.0 * dz_z(x[0], x[1], pot)) * x[6] - std::exp(2 * (log_k - x[2])) * x[4];
}

void _equations_tensor(double dx_dt[], const double, const double x[], void* data)
{
    auto ptr = static_cast<void**>(data);
    auto pot = static_cast<Potential*> (ptr[0]);
    auto params = static_cast<double*> (ptr[1]);
    double log_k = std::log(params[2]);
    
    dx_dt[0] = x[1];
    dx_dt[1] = - (3 * _H(x[0], x[1], pot) * x[1] + pot->dV(x[0]));
    dx_dt[2] = _H(x[0], x[1], pot);
    dx_dt[3] = x[5];
    dx_dt[4] = x[6];
    dx_dt[5] = -(3 * _H(x[0], x[1], pot)) * x[5] - std::exp(2 * (log_k - x[2])) * x[3];
    dx_dt[6] = -(3 * _H(x[0], x[1], pot)) * x[6] - std::exp(2 * (log_k - x[2])) * x[4];
}

void _inflation_end(double g[], const double, const double x[], void* data)
{
    auto ptr = static_cast<void**>(data);
    auto pot = static_cast<Potential*> (ptr[0]);
    auto params = static_cast<double*> (ptr[1]);
    
    auto p = -(2.5 * x[1] * x[1] + x[1] * pot->dV(x[0]) / _H(x[0], x[1], pot) - pow(0.5 * x[1] * x[1] / _H(x[0], x[1], pot), 2));
    if(p >= 0)
        g[0] = (-_H(x[0], x[1], pot) + 0.5 * x[1] * x[1] / _H(x[0], x[1], pot));
    else if(p < 0)
        g[0] = p;
}

void _inflation_begin(double g[], const double, const double x[], void* data)
{
    auto ptr = static_cast<void**>(data);
    auto pot = static_cast<Potential*> (ptr[0]);
    
    g[0] = (-_H(x[0], x[1], pot) + 0.5 * x[1] * x[1] / _H(x[0], x[1], pot));
    
}

void _Find_N(double g[], const double, const double x[], void* data)
{
    auto ptr = static_cast<void**>(data);
    auto params = static_cast<double*> (ptr[1]);
    
    g[0] = params[0] - x[2];
}

void _start(double g[], const double, const double x[], void* data)
{
    auto ptr = static_cast<void**>(data);
    auto pot = static_cast<Potential*> (ptr[0]);
    auto params = static_cast<double*> (ptr[1]);
    
    g[0] = (params[0] + std::log(200.)) + (x[2] + std::log(_H(x[0], x[1], pot)));
}

void _finish(double g[], const double, const double x[], void* data)
{
    auto ptr = static_cast<void**>(data);
    auto pot = static_cast<Potential*> (ptr[0]);
    auto params = static_cast<double*> (ptr[1]);
    
    g[0] = (params[0] - std::log(200.)) + (x[2] + std::log(_H(x[0], x[1], pot)));
}

NumericModeSolver::NumericModeSolver(Potential* _pot, double _N_star, double _N_r) : 
    pot{_pot}, N_star{_N_star}, N_dagger{}, N_r{_N_r}, phi_p{}, dphi_p{}, log_aH_star{}, phi_IC{}, dphi_IC{}, n_IC{}, N_end{}
{
    void* ptrs[2];
    ptrs[0] = static_cast<void*> (pot);
    double params[3];
    params[0] = params[1] = params[2] = 0;
    ptrs[1] = static_cast<void*> (params);
    
    double t0 = 1;
    double N_temp = 0;
    phi_p = 0;
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
    
    dphi_p = - pot->dV(phi_p) / (3 * std::sqrt(pot->V(phi_p) / 3));
    
    //Find N_end
    double t = t0;
    std::vector<double> x = {phi_p, dphi_p, 0., 0.};
    dlsodar desolver(4, 1, 1000000);
    desolver.integrate(t, 1e10, &x[0], __equations, _inflation_end, static_cast<void*> (ptrs));
    N_end = x[2];
    params[0] = N_end;
    
    //Find aH_star
    params[0] = N_end - N_star;
    t = t0;
    x = {phi_p, dphi_p, 0., 0.};
    desolver = dlsodar(4, 1, 1000000);
    desolver.integrate(t, 1e10, &x[0], __equations, _Find_N, static_cast<void*> (ptrs));
    log_aH_star = x[2] + std::log(_H(x[0], x[1], pot));
    
    //Find BackgroundVar at N_set
    params[0] = N_end - N_r;
    t = t0;
    x = {phi_p, dphi_p, 0., 0.};
    desolver = dlsodar(4, 1, 1000000);
    desolver.integrate(t, 1e10, &x[0], __equations, _Find_N, static_cast<void*> (ptrs));
    phi_IC = x[0];
    dphi_IC = x[1];
    n_IC = x[2];
    
}

NumericModeSolver::NumericModeSolver(Potential* _pot, double _N_star, double _N_dagger, double _N_r) :
    pot{_pot}, N_star{_N_star}, N_dagger{_N_dagger}, N_r{_N_r}, phi_p{}, dphi_p{}, log_aH_star{}, phi_IC{}, dphi_IC{}, n_IC{}, N_end{}
{
    void* ptrs[2];
    ptrs[0] = static_cast<void*> (pot);
    double params[3];
    params[0] = params[1] = params[2] = 0;
    ptrs[1] = static_cast<void*> (params);
    
    double t0 = 1;
    phi_p = 0, dphi_p = - std::sqrt(2.0/3), N_end = 0;
    double N_temp_d = 0, NN = 0;
    while(N_end < N_star + N_dagger)
    {
        phi_p += 0.5;
        double t = t0;
        std::vector<double> x = {phi_p, dphi_p, 0. ,0.};
        dlsodar desolver(4, 1, 1000000);
        desolver.integrate(t, 1e10, &x[0], __equations, _inflation_end, static_cast<void*> (ptrs));
        N_end = x[2];
    }
    double phi_old = phi_p - 0.5;
    phi_p += 0.001;
    while(abs(N_temp_d - N_dagger) > 0.01)
    {
        double t = t0;
        std::vector<double> x = {phi_p, dphi_p, 0., 0.};
        dlsodar desolver(4, 1, 1000000);
        desolver.integrate(t, 1e10, &x[0], __equations, _inflation_end, static_cast<void*> (ptrs));
        N_end = x[2];
        
        t = t0;
        x = {phi_p, dphi_p, 0., 0.};
        desolver = dlsodar(4, 1, 100000);
        desolver.integrate(t, 1e10, &x[0], __equations, _inflation_begin, static_cast<void*> (ptrs));
        NN = x[2];
        
        params[0] = N_end - N_star;
        ptrs[1] = static_cast<void*> (params);
        desolver.integrate(t, 1e10, &x[0], __equations, _Find_N, static_cast<void*> (ptrs));
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
    
    //Find aH_star
    params[0] = N_end - N_star;
    auto t = t0;
    std::vector<double> x = {phi_p, dphi_p, 0., 0.};
    dlsodar desolver(4, 1, 1000000);
    desolver.integrate(t, 1e10, &x[0], __equations, _Find_N, static_cast<void*> (ptrs));
    log_aH_star = x[2] + std::log(_H(x[0], x[1], pot));
    
    //Find BackgroundVar at N_set
    params[0] = N_end - N_r;
    t = t0;
    x = {phi_p, dphi_p, 0., 0.}; 
    desolver = dlsodar(4, 1, 1000000);
    desolver.integrate(t, 1e10, &x[0], __equations, _Find_N, static_cast<void*> (ptrs));
    phi_IC = x[0];
    dphi_IC = x[1];
    n_IC = x[2];
}


double NumericModeSolver::Find_PPS_Scalar(double k)
{
    void* ptrs[2];
    ptrs[0] = static_cast<void*> (pot);
    double params[3];
    params[0] = params[1] = params[2] = 0;
    ptrs[1] = static_cast<void*> (params);
    
    double t0 = 1;
    //Find R1, R2 at eta_r for matching
    double R1_IC, R2_IC, dR1_IC, dR2_IC;
    //Find R1, R2 at end of inflation
    double R1, R2, dR1, dR2;
    
    //Solve MS
    k *= std::exp(log_aH_star) / (0.05);
    params[2] = k;
    double t = t0;
    std::vector<double> x = {phi_p, dphi_p, 0, 1, 0, 0, 1};
    
    dlsodar desolver(7, 1, 100000000);
    params[0] = N_end - N_r;
    desolver.integrate(t, 1e10, &x[0], _equations, _Find_N, static_cast<void*> (ptrs));
    R1_IC = x[3];
    R2_IC = x[4];
    dR1_IC = x[5];
    dR2_IC = x[6];
    
    desolver.integrate(t, 1e10, &x[0], _equations, _inflation_end, static_cast<void*> (ptrs));
    R1 = x[3];
    R2 = x[4];
    dR1 = x[5];
    dR2 = x[6];
    
    //Matching and find PPS
    double PPS;
    
    double Z = std::exp(n_IC) * dphi_IC / _H(phi_IC, dphi_IC, pot);
    double dZ_Z = dz_z(phi_IC, dphi_IC, pot);
    
    //Set Vacuum Initial Conditions
    auto R0 = 1.0 / (Z * std::sqrt(2 * k));
    auto dR0 = 0;//(-1.0 * I * k * std::exp(- n_IC) - dZ_Z) * R0;
    
    auto det = R1_IC * dR2_IC - R2_IC * dR1_IC;
    auto A = (R0 * dR2_IC - R2_IC * dR0) / det;
    auto B = (-R0 * dR1_IC + R1_IC * dR0) / det;
    
    PPS = (pow(k, 3) / (2 * M_PI * M_PI)) * pow(abs(A * R1 + B * R2), 2);
    
    return PPS;
}

double NumericModeSolver::Find_PPS_Tensor(double k)
{
    void* ptrs[2];
    ptrs[0] = static_cast<void*> (pot);
    double params[3];
    params[0] = params[1] = params[2] = 0;
    ptrs[1] = static_cast<void*> (params);
    
    double t0 = 1;
    //Find R1, R2 at eta_r for matching
    double R1_IC, R2_IC, dR1_IC, dR2_IC;
    //Find R1, R2 at end of inflation
    double R1, R2, dR1, dR2;
    
    //Solve MS
    k *= std::exp(log_aH_star) / (0.05);
    params[2] = k;
    double t = t0;
    std::vector<double> x = {phi_p, dphi_p, 0, 1, 0, 0, 1};
    
    dlsodar desolver(7, 1, 100000000);
    params[0] = N_end - N_r;
    desolver.integrate(t, 1e10, &x[0], _equations_tensor, _Find_N, static_cast<void*> (ptrs));
    R1_IC = x[3];
    R2_IC = x[4];
    dR1_IC = x[5];
    dR2_IC = x[6];
    
    desolver.integrate(t, 1e10, &x[0], _equations_tensor, _inflation_end, static_cast<void*> (ptrs));
    R1 = x[3];
    R2 = x[4];
    dR1 = x[5];
    dR2 = x[6];
    
    //Matching and find PPS
    double PPS;
    
    //Set Vacuum Initial Conditions
    auto R0 = 2.0 / (std::exp(n_IC) * std::sqrt(2 * k));
    auto dR0 = (-1.0 * I * k * std::exp(- n_IC) - _H(phi_IC, dphi_IC, pot)) * R0;
    
    auto det = R1_IC * dR2_IC - R2_IC * dR1_IC;
    auto A = (R0 * dR2_IC - R2_IC * dR0) / det;
    auto B = (-R0 * dR1_IC + R1_IC * dR0) / det;
    
    PPS = (pow(k, 3) / (2 * M_PI * M_PI)) * pow(abs(A * R1 + B * R2), 2);
    
    return PPS;
}

double NumericModeSolver::Find_PPS(double k)
{
    void* ptrs[2];
    ptrs[0] = static_cast<void*> (pot);
    double params[3];
    params[0] = params[1] = params[2] = 0;
    ptrs[1] = static_cast<void*> (params);
    
    double t0 = 1;
    double N_temp = 0;
    phi_p = 0;
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
    
    dphi_p = - pot->dV(phi_p) / (3 * std::sqrt(pot->V(phi_p) / 3));
    
    //Find N_end
    double t = t0;
    std::vector<double> x = {phi_p, dphi_p, 0., 0.};
    dlsodar desolver(4, 1, 1000000);
    desolver.integrate(t, 1e10, &x[0], __equations, _inflation_end, static_cast<void*> (ptrs));
    N_end = x[2];
    params[0] = N_end;
    
    //Find aH_star
    params[0] = N_end - N_star;
    t = t0;
    x = {phi_p, dphi_p, 0., 0.};
    desolver = dlsodar(4, 1, 1000000);
    desolver.integrate(t, 1e10, &x[0], __equations, _Find_N, static_cast<void*> (ptrs));
    log_aH_star = x[2] + std::log(_H(x[0], x[1], pot));
    
    k *= std::exp(log_aH_star) / (0.05);
    
    //Find BackgroundVar at N_start
    params[0] = -std::log(k);
    t = t0;
    x = {phi_p, dphi_p, 0., 0.};
    desolver = dlsodar(4, 1, 1000000);
    desolver.integrate(t, 1e10, &x[0], __equations, _start, static_cast<void*> (ptrs));
    auto t_start = t;
    auto phi_start = x[0];
    auto dphi_start = x[1];
    auto n_start = x[2];
    auto eta_start = x[3];
    
    //Find BackgroundVar at N_finish
    params[0] = -std::log(k);
    x = {phi_p, dphi_p, 0., 0.};
    desolver = dlsodar(4, 1, 1000000);
    desolver.integrate(t, 1e10, &x[0], __equations, _finish, static_cast<void*> (ptrs));
    auto t_finish = t;
    auto phi_finish = x[0];
    auto dphi_finish = x[1];
    auto n_finish = x[2];
    auto eta_finish = x[3];
    
    //Find R1, R2 at eta_r for matching
    double R1_IC = 1, R2_IC = 0, dR1_IC = 0, dR2_IC = 1;
    //Find R1, R2 at end of inflation
    double R1, R2, dR1, dR2;
    
    //Solve MS
    params[2] = k;
    x = {phi_start, dphi_start, n_start, 1, 0, 0, 1};
    desolver = dlsodar(7, 1, 1000000);
    desolver.integrate(t_start, t_finish, &x[0], _equations, _inflation_end, static_cast<void*> (ptrs));
    R1 = x[3];
    R2 = x[4];
    dR1 = x[5];
    dR2 = x[6];
    
    //Matching and find PPS
    double PPS;
    double Z = std::exp(n_start) * dphi_start / _H(phi_start, dphi_start, pot);
    double dZ_Z = dz_z(phi_start, dphi_start, pot);
    
    //Set Vacuum Initial Conditions
    auto R0 = (cos(k*eta_start) + I * sin(k*eta_start)) / (Z * std::sqrt(2 * k));
    auto dR0 = (-1.0 * I * k * std::exp(-n_start) - dZ_Z) * R0;
    
    PPS = (pow(k, 3) / (2 * M_PI * M_PI)) * pow(abs(R0 * R1 + dR0 * R2), 2);
    
    return PPS;
}
