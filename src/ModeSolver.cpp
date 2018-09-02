#include "ModeSolver.hpp"

ModeSolver::ModeSolver(BackgroundSolution _Bsol): BasicModeSolver(),
Bsol{_Bsol}, Tsol{}, N_r{0}, PPS_error{5e-3}, vacuum{BD} {}

void ModeSolver::Initial_Conditions(VacuumChoice _vacuum, double _N_r)
{
    vacuum = _vacuum;
    N_r = Bsol.N_end - _N_r;
}

double ModeSolver::Find_PPS_Scalar(double k)
{
    //Scale k
    k *= Bsol.aH_star / 0.05;
    
    //Find Transitions
    double N_f = log(k / 0.05) + 55;
    if(N_f > Bsol.N_end)
    {
        N_f = Bsol.N_end;
    }
    Transitions T(N_r, N_f, Bsol.omega_2, Bsol.log_aH, Bsol.N_extrema);
    Tsol = T.Find(k, exp((log(PPS_error) - 1.0633) / 0.7852));
    
    //Match and evaluate PPS
    Eigen::Vector2cd Q = Match(k);
    double F = exp(N_f) * Bsol.dphi_H(N_f) * exp(0.5 * Bsol.log_aH(N_f));
    
    return (pow(k, 3) / (2 * M_PI * M_PI)) * pow(abs(Q[0] / F), 2);
}

double ModeSolver::Find_PPS_Tensor(double k)
{
    //Scale k
    k *= Bsol.aH_star / 0.05;
    
    //Find Transitions
    double N_f = log(k / 0.05) + 55;
    if(N_f > Bsol.N_end)
    {
        N_f = Bsol.N_end;
    }
    Transitions T(N_r, N_f, Bsol.omega_2_tensor, Bsol.log_aH, Bsol.N_extrema_tensor);
    Tsol = T.Find(k, exp((log(PPS_error) + 0.4179) / 0.5615));

    //Match and evaluate PPS
    Eigen::Vector2cd Q = Match(k);
    double F = exp(N_f) * exp(0.5 * Bsol.log_aH(N_f));
    
    return 4 * (pow(k, 3) / (2 * M_PI * M_PI)) * pow(abs(Q[0] / F), 2);
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
    
    //Evolution matrix
    Eigen::Matrix2cd Evolve = Eigen::Matrix2d::Identity();
    for(size_t n = 0; n < Tsol.seg_control.size()-1; n++)
    {
        if(Tsol.seg_control[n].second == 0)
        {
            Evolve = Airy_Mat(Tsol.a[n], Tsol.b[n], Tsol.seg_control[n].first, Tsol.seg_control[n+1].first) * Evolve;
        }
        else if(Tsol.seg_control[n].second == 1)
        {
            Evolve = Bessel_Mat(Tsol.a[n], Tsol.b[n], Tsol.seg_control[n].first, Tsol.seg_control[n+1].first) * Evolve;
        }
        else if(Tsol.seg_control[n].second == -1)
        {
            Evolve = Modified_Bessel_Mat(Tsol.a[n], Tsol.b[n], Tsol.seg_control[n].first, Tsol.seg_control[n+1].first) * Evolve;
        }
    }
    
    return Evolve * Q;
}

Eigen::Matrix2d ModeSolver::Airy_Mat(double a, double b, double N0, double N1)
{
    Eigen::Matrix2d A;

    double p = pow(abs(b), 1.0/3.0);
    double x0 = -((a + b * N0) /p/p);
    double x1 = -((a + b * N1) /p/p);
    
    if(b < 0)
    {
        A = Airy_gen(p, x0, x1);
    }
    else if(b > 0)
    {
        A = Airy_gen(-p, x0, x1);
    }
    
    return A;
}

Eigen::Matrix2d ModeSolver::Bessel_Mat(double a, double b, double N0, double N1)
{
    double p = b;
    double x0 = exp(a + b * N0);
    double x1 = exp(a + b * N1);
    
    Eigen::Matrix2d B = Bessel_gen(p, x0, x1);

    return B;
}

Eigen::Matrix2d ModeSolver::Modified_Bessel_Mat(double a, double b, double N0, double N1)
{
    double p = b;
    double x0 = exp(a + b * N0);
    double x1 = exp(a + b * N1);
    
    Eigen::Matrix2d MB = Modified_Bessel_gen(p, x0, x1);
    
    return MB;
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

Eigen::Matrix2d ModeSolver::Bessel_gen(double p, double x0, double x1)
{
    Eigen::Matrix2cd B0, B1;
    Eigen::Matrix2d B;
    
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
    
    B = (B1 * B0.inverse()).real();
    
    return B;
}

Eigen::Matrix2d ModeSolver::Modified_Bessel_gen(double p, double x0, double x1)
{
    Eigen::Matrix2cd MB0, MB1;
    Eigen::Matrix2d MB;
    
    if((x0 / p) < 20 and (x1 / p) < 20 and (x0 / p) > -15 and (x1 / p) > -15)
    {
        auto I0 = Bessel_I(0, x0 / p);
        auto K0 = Bessel_K(0, x0 / p);
        auto Ip0 = x0 * Bessel_I(1, x0 / p);
        auto Kp0 = -x0 * Bessel_K(1, x0 / p);
        
        auto I = Bessel_I(0, x1 / p);
        auto K = Bessel_K(0, x1 / p);
        auto Ip = x1 * Bessel_I(1, x1 / p);
        auto Kp = -x1 * Bessel_K(1, x1 / p);
        
        MB0 << I0,  K0,
            Ip0, Kp0;
        
        MB1 << I,  K,
            Ip, Kp;

        MB = (MB1 * MB0.inverse()).real();
    }
    else
    {
        auto a00 = 0.5 * (exp((x1-x0) / p) + exp(-(x1-x0) / p)) * pow(x0/x1, 0.5);
        auto a01 = 0.5 * (exp((x1-x0) / p) - exp(-(x1-x0) / p)) * pow(x0 * x1, -0.5);
        auto a10 = 0.5 * (exp((x1-x0) / p) - exp(-(x1-x0) / p)) * pow(x0 * x1, 0.5);
        auto a11 = 0.5 * (exp((x1-x0) / p) + exp(-(x1-x0) / p)) * pow(x0/x1, -0.5);
        
        MB<< a00, a01,
            a10, a11;
    }
    
    return MB;
}


