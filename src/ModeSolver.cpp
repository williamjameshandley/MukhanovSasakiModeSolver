#include "ModeSolver.hpp"

ModeSolver::ModeSolver(BackgroundSolution _Bsol):
Bsol{_Bsol}, Tsol{}, N_r{0}, vacuum{BD}, PPS{} {}

void ModeSolver::Initial_Conditions(VacuumChoice _vacuum, double _N_r)
{
    vacuum = _vacuum;
    N_r = _N_r;
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

double ModeSolver::Find_PPS(double k)
{
    k_plot.push_back(k);
    k *= Bsol.aH_star / 0.05;
    
    double N_f = log(k / 0.05) + 55;
    Transitions T(N_r, N_f, Bsol);
    Tsol = T.Find(k, 1e-4);
    
    Eigen::Vector2cd Q = Match(k);
    
    double F = exp(N_f) * Bsol.dphi_H(N_f) * exp(0.5 * Bsol.log_aH(N_f));
    
    return (pow(k, 3) / (2 * M_PI * M_PI)) * pow(abs(Q[0] / F), 2);
}

void ModeSolver::Construct_PPS(double k_i, double k_f, double error = 1e-3)
{
    std::vector<std::pair<double, double>> k_pair;
    
    k_pair.push_back(std::make_pair(k_i, k_f));
    
    PPS.insert(k_pair[0].first, Find_PPS(k_pair[0].first));
    PPS.insert(k_pair[0].second, Find_PPS(k_pair[0].second));
    
    double lim = error;
    while(k_pair.size() != 0)
    {
        for(size_t n = 0; n < k_pair.size(); n++)
        {
            k_i = k_pair[n].first;
            k_f = k_pair[n].second;
            
            auto k_m1 = exp((2 * log(k_i) + log(k_f)) / 3.0);
            auto k_m2 = exp((log(k_i) + 2 * log(k_f)) / 3.0);
            
            auto temp_approx1 = (PPS(k_m1));
            auto temp_approx2 = (PPS(k_m2));
            auto temp_true1 = (Find_PPS(k_m1));
            auto temp_true2 = (Find_PPS(k_m2));
            
            k_pair.erase(k_pair.begin() + static_cast<int>(n));
            
            if(abs(temp_true1 - temp_approx1) / abs(temp_true1) > lim and abs(temp_true2 - temp_approx2) / abs(temp_true2) > lim)
            {
                k_pair.insert(k_pair.begin() + static_cast<int>(n), std::make_pair(k_i, k_m1));
                k_pair.insert(k_pair.begin() + static_cast<int>(n) + 1, std::make_pair(k_m1, k_m2));
                k_pair.insert(k_pair.begin() + static_cast<int>(n) + 2, std::make_pair(k_m2, k_f));
                n += 2;
                PPS.insert(k_m1, (temp_true1));
                PPS.insert(k_m2, (temp_true2));
                k_plot.push_back(k_m1);
                k_plot.push_back(k_m2);
            }
            if(abs(temp_true1 - temp_approx1) / abs(temp_true1) > lim and abs(temp_true2 - temp_approx2) / abs(temp_true2) < lim)
            {
                k_pair.insert(k_pair.begin() + static_cast<int>(n), std::make_pair(k_i, k_m1));
                k_pair.insert(k_pair.begin() + static_cast<int>(n) + 1, std::make_pair(k_m1, k_m2));
                n += 1;
                PPS.insert(k_m1, (temp_true1));
                k_plot.push_back(k_m1);
            }
            if(abs(temp_true1 - temp_approx1) / abs(temp_true1) < lim and abs(temp_true2 - temp_approx2) / abs(temp_true2) > lim)
            {
                k_pair.insert(k_pair.begin() + static_cast<int>(n), std::make_pair(k_m1, k_m2));
                k_pair.insert(k_pair.begin() + static_cast<int>(n) + 1, std::make_pair(k_m2, k_f));
                n += 1;
                PPS.insert(k_m2, (temp_true2));
                k_plot.push_back(k_m2);
            }
        }
    }
    std::sort(k_plot.begin(), k_plot.end());
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
