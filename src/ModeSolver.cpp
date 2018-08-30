#include "ModeSolver.hpp"

ModeSolver::ModeSolver(BackgroundSolution _Bsol): BasicModeSolver(),
Bsol{_Bsol}, Tsol{}, N_r{0}, vacuum{BD} {}

void ModeSolver::Initial_Conditions(VacuumChoice _vacuum, double _N_r)
{
    vacuum = _vacuum;
    N_r = Bsol.N_end - _N_r;
}

Eigen::Vector2cd ModeSolver::Match(double k)
{
    //Setting Vacuum
    double aH_r = exp(Bsol.log_aH(N_r));
    double epsilon = 0.5 * pow(Bsol.dphi_H(N_r), 2);
    Eigen::Vector2cd Q;

    Q[0] = sqrt(aH_r / (2 * k));
    
    if(vacuum == BD)
        Q[1] = Q[0] * (0.5 * (1 - epsilon) + (-I * k / aH_r));
    else
        throw std::runtime_error("Initial conditions unknown");
    
    //Match
    return Airy_Mat() * Bessel_Mat() * Q;
}

double ModeSolver::Find_PPS_Scalar(double k)
{
    k *= Bsol.aH_star / 0.05;
    
    double N_f = log(k / 0.05) + 55;
    if(N_f > Bsol.N_end)
    {
        N_f = Bsol.N_end;
    }
    Transitions T(N_r, N_f, Bsol.omega_2, Bsol.log_aH, Bsol.N_extrema);
    Tsol = T.Find(k, 1e-4);
    
    Eigen::Vector2cd Q = Match(k);
    
    double F = exp(N_f) * Bsol.dphi_H(N_f) * exp(0.5 * Bsol.log_aH(N_f));
    
    return (pow(k, 3) / (2 * M_PI * M_PI)) * pow(abs(Q[0] / F), 2);
}

double ModeSolver::Find_PPS_Tensor(double k)
{
    k *= Bsol.aH_star / 0.05;
    
    double N_f = log(k / 0.05) + 55;
    if(N_f > Bsol.N_end)
    {
        N_f = Bsol.N_end;
    }
    Transitions T(N_r, N_f, Bsol.omega_2_tensor, Bsol.log_aH, Bsol.N_extrema_tensor);
    Tsol = T.Find(k, 1e-4);
    
    Eigen::Vector2cd Q = Match(k);
    
    double F = exp(N_f) * exp(0.5 * Bsol.log_aH(N_f));
    
    return 4 * (pow(k, 3) / (2 * M_PI * M_PI)) * pow(abs(Q[0] / F), 2);
}

Eigen::Matrix2d ModeSolver::Airy_Mat()
{
    Eigen::Matrix2d Mat = Eigen::Matrix2d::Identity();
    
    for(size_t n = 0; n < Tsol.lin_a.size(); n++)
    {
        double p = pow(abs(Tsol.lin_b[n]), 1.0/3.0);
        double x0 = -((Tsol.lin_a[n] + Tsol.lin_b[n] * Tsol.lin_N_step[n]) /p/p);
        double x1 = -((Tsol.lin_a[n] + Tsol.lin_b[n] * Tsol.lin_N_step[n+1]) /p/p);
        
        if(Tsol.lin_b[n] < 0)
        {
            Eigen::Matrix2d A = Airy_gen(p, x0, x1);
            Mat = A * Mat;
        }
        else if(Tsol.lin_b[n] > 0)
        {
            Eigen::Matrix2d A = Airy_gen(-p, x0, x1);
            Mat = A * Mat;
        }
    }
    
    return Mat;
}

Eigen::Matrix2cd ModeSolver::Bessel_Mat()
{
    Eigen::Matrix2cd Mat = Eigen::Matrix2d::Identity();
    
    if(Tsol.log_a.size() != 0)
    {
        for(size_t n = 0; n < Tsol.log_a.size(); n++)
        {
            double p = Tsol.log_b[n];
            double x0 = exp(Tsol.log_a[n] + Tsol.log_b[n] * Tsol.log_N_step[n]);
            double x1 = exp(Tsol.log_a[n] + Tsol.log_b[n] * Tsol.log_N_step[n+1]);
            
            Eigen::Matrix2cd B = Bessel_gen(p, x0, x1);
            
            Mat = B * Mat;
        }
    }
    return Mat;
}

Eigen::Matrix2d ModeSolver::Airy_gen(double p, double x0, double x1)
{
    Eigen::Matrix2d A0, A1, A;
    double Ai0, Bi0, Aip0, Bip0, Ai, Bi, Aip, Bip;
    
    if(x0 < 25 and x1 < 25)
    {
        Airy(x0, Ai0, Aip0, Bi0, Bip0);
        Airy(x1, Ai, Aip, Bi, Bip);
        
        Aip0 *= p;
        Bip0 *= p;
        Aip *= p;
        Bip *= p;
        
        A0 << Ai0,  Bi0,
              Aip0, Bip0;
        
        A1 << Ai,  Bi,
              Aip, Bip;
        
        A = A1 * A0.inverse();
    }
    else
    {
        auto a0 = pow(x0, 1.5), b0 = pow(x0, 0.25), a1 = pow(x1, 1.5), b1 = pow(x1, 0.25);
        
        auto a00 = 0.5 * (exp((2.0/3) * (a1 - a0)) + exp(-(2.0/3) * (a1 - a0))) * b0 / b1;
        auto a01 = 0.5 * (exp((2.0/3) * (a1 - a0)) - exp(-(2.0/3) * (a1 - a0))) / (p * b0 * b1);
        auto a10 = 0.5 * (exp((2.0/3) * (a1 - a0)) - exp(-(2.0/3) * (a1 - a0))) * (p * b0 * b1);
        auto a11 = 0.5 * (exp((2.0/3) * (a1 - a0)) + exp(-(2.0/3) * (a1 - a0))) * b1 / b0;
        
        A << a00,  a01,
             a10, a11;
    }
    
    return A;
}

Eigen::Matrix2cd ModeSolver::Bessel_gen(double p, double x0, double x1)
{
    Eigen::Matrix2cd B0, B1;
    
    auto J0 = Bessel_J(0, x0 / p);
    auto Y0 = Bessel_Y(0, x0 / p);
    auto Jp0 = -x0 * Bessel_J(1, x0 / p);
    auto Yp0 = -x0 * Bessel_Y(1, x0 / p);
    
    auto J = Bessel_J(0, x1 / p);
    auto Y = Bessel_Y(0, x1 / p);
    auto Jp = -x1 * Bessel_J(1, x1 / p);
    auto Yp = -x1 * Bessel_Y(1, x1 / p);
    
    B0 << J0,  Y0,
          Jp0, Yp0;
    
    B1 << J,  Y,
          Jp, Yp;
    
    return B1 * B0.inverse();
}

