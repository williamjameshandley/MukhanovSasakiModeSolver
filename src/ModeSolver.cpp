#include "ModeSolver.hpp"

ModeSolver::ModeSolver(BackgroundSolution _Bsol):
Bsol{_Bsol}, Tsol{}, N_r{0}, vacuum{BD}, Z{}, H{}, DPHI{}, PPS{}
{
    for(size_t o = 0; o < Bsol.N.size(); o++)
    {
        DPHI.insert(Bsol.N[o], Bsol.dphi[o]);
        Z.insert(Bsol.N[o], Bsol.z[o]);
        H.insert(Bsol.N[o], Bsol.H[o]);
    }
}

void ModeSolver::Initial_Conditions(VacuumChoice _vacuum, double _N_r)
{
    vacuum = _vacuum;
    N_r = _N_r;
}

Eigen::Matrix2d ModeSolver::Mat()
{
    Eigen::Matrix2d Mat = Eigen::Matrix2d::Identity();
    
    for(size_t n = Tsol.N_step.size() - 2; n != static_cast<size_t>(-1); n--)
    {
        double p = pow(abs(Tsol.b[n]), 1.0/3.0);
        double x0 = -((Tsol.a[n] + Tsol.b[n] * Tsol.N_step[n]) /p/p);
        double x1 = -((Tsol.a[n] + Tsol.b[n] * Tsol.N_step[n+1]) /p/p);
        
        Eigen::Matrix2d A = Airy_mat(p, x0, x1);
        
        Mat *= A;
    }

    return Mat;
}

Eigen::Vector2cd ModeSolver::Match(double k)
{
    //Setting Vacuum
    double aH_r = exp(N_r) * H(N_r);
    double dphi_H_2 = pow(0.5 * DPHI(N_r) / H(N_r), 2);
    Eigen::Vector2cd Q;
    
    Q[0] = sqrt(aH_r / (2 * k));
    
    if(vacuum == BD)
        Q[1] = Q[0] * (0.5 - dphi_H_2 + (-I * k / aH_r));
        //only change the last part for different IC (only the -ik part)
    else
        throw std::runtime_error("Initial conditions unknown");

    //Match
    return Mat() * Q;
}

double ModeSolver::Find_PPS(double k)
{
    double N_f = 0.98 * Bsol.N.back();
    Transitions T(N_r, N_f, Bsol);
    Tsol = T.Find(k, 1e-4);
    
    Eigen::Vector2cd Q = Match(k);
    
    double F = Z(N_f) * sqrt(exp(N_f) * H(N_f));
    
    return (pow(k, 3) / (2 * M_PI * M_PI)) * pow(abs(Q[0] / F), 2);
}

void ModeSolver::Construct_PPS(double k0, double k1)
{
    std::vector<std::pair<double, double>> k_pair;
    
    k_pair.push_back(std::make_pair(k0, k1));
    
    PPS.insert(k_pair[0].first, Find_PPS(k_pair[0].first));
    PPS.insert(k_pair[0].second, Find_PPS(k_pair[0].second));
    
    double lim = 1e-3;
    
    while(k_pair.size() != 0)
    {
        for(size_t n = 0; n < k_pair.size(); n++)
        {
            k0 = k_pair[n].first;
            k1 = k_pair[n].second;
            
            auto k_m1 = exp((2 * log(k0) + log(k1)) / 3.0);
            auto k_m2 = exp((log(k0) + 2 * log(k1)) / 3.0);
            
            auto temp_approx1 = (PPS(k_m1));
            auto temp_approx2 = (PPS(k_m2));
            auto temp_true1 = (Find_PPS(k_m1));
            auto temp_true2 = (Find_PPS(k_m2));

            k_pair.erase(k_pair.begin() + static_cast<int>(n));
            
            if(abs(temp_true1 - temp_approx1) / temp_true1 > lim and abs(temp_true2 - temp_approx2) / temp_true2 > lim)
            {
                k_pair.insert(k_pair.begin() + static_cast<int>(n), std::make_pair(k0, k_m1));
                k_pair.insert(k_pair.begin() + static_cast<int>(n) + 1, std::make_pair(k_m1, k_m2));
                k_pair.insert(k_pair.begin() + static_cast<int>(n) + 2, std::make_pair(k_m2, k1));
                n += 2;
            }
            else if(abs(temp_true1 - temp_approx1) / temp_true1 > lim and abs(temp_true2 - temp_approx2) / temp_true2 < lim)
            {
                k_pair.insert(k_pair.begin() + static_cast<int>(n), std::make_pair(k0, k_m1));
                k_pair.insert(k_pair.begin() + static_cast<int>(n) + 1, std::make_pair(k_m1, k1));
                n += 1;
            }
            else if(abs(temp_true1 - temp_approx1) / temp_true1 < lim and abs(temp_true2 - temp_approx2) / temp_true2 > lim)
            {
                k_pair.insert(k_pair.begin() + static_cast<int>(n), std::make_pair(k0, k_m2));
                k_pair.insert(k_pair.begin() + static_cast<int>(n) + 1, std::make_pair(k_m2, k1));
                n += 1;
            }
            PPS.insert(k_m1, (temp_true1));
            PPS.insert(k_m2, (temp_true2));
        }
    }
}

Eigen::Matrix2d ModeSolver::Airy_mat(double p, double x0, double x1)
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

