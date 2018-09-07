#include "ModeSolver.hpp"

ModeSolver::ModeSolver(BackgroundSolution _Bsol): Bsol{_Bsol}, N_r{0}, PPS_error{5e-3}, vacuum{BD} {}

void ModeSolver::Initial_Conditions(VacuumChoice _vacuum, double _N_r)
{
    vacuum = _vacuum;
    N_r = Bsol.N_end - _N_r;
}

double ModeSolver::Find_PPS_Scalar(double k)
{
    //Scale k
    k *= Bsol.aH_star / 0.05;
    
    //Initial Q
    Eigen::Vector2cd Q_i = Initial_Q(k);
    
    //Find N_f
    double N_f = log(k / 0.05) + 55;
    if(N_f > Bsol.N_end)
        N_f = Bsol.N_end;
    
    //Evolve Q
    Eigen::Vector2cd Q_f = Evolve(Q_i, k, N_r, N_f);
    
    double F = exp(N_f) * Bsol.dphi_H(N_f) * exp(0.5 * Bsol.log_aH(N_f));
    
    return (pow(k, 3) / (2 * M_PI * M_PI)) * pow(abs(Q_f[0] / F), 2);
}

double ModeSolver::Find_PPS_Tensor(double k)
{
    return 0;
}

Eigen::Vector2cd ModeSolver::Initial_Q(double k)
{
    //Setting Vacuum
    double aH_r = exp(Bsol.log_aH(N_r));
    double epsilon = 0.5 * pow(Bsol.dphi_H(N_r), 2);
    Eigen::Vector2cd Q_i;
    
    Q_i[0] = sqrt(aH_r / (2 * k));
    
    if(vacuum == BD)
        Q_i[1] = Q_i[0] * (0.5 * (1 - epsilon) + (-I * k / aH_r));
    else
        throw std::runtime_error("Initial conditions unknown");
    
    return Q_i;
}

Eigen::Vector2cd ModeSolver::Evolve(Eigen::Vector2cd Q_i, double k, double N_initial, double N_final)
{
    std::map<double, Vars> Seg;
    Vars Null_var(0, 0, 0, Eigen::Matrix2d::Identity());
    
    //Initialize Extrema and end points
    Seg.insert(std::pair<double,Vars> (N_initial, Null_var));
    Seg.insert(std::pair<double,Vars> (N_final, Null_var));
    if(Bsol.N_extrema.size() != 0)
    {
        auto p0 = static_cast<size_t>(std::lower_bound(Bsol.N_extrema.begin(), Bsol.N_extrema.end(), N_initial) - Bsol.N_extrema.begin());
        auto p1 = static_cast<size_t>(std::lower_bound(Bsol.N_extrema.begin(), Bsol.N_extrema.end(), N_final) - Bsol.N_extrema.begin());
        if(p1 <= Bsol.N_extrema.size() and p1 > p0)
        {
            for(size_t n = p0; n < p1; n++)
                Seg.insert(std::pair<double,Vars> (Bsol.N_extrema[n], Null_var));
        }
    }
    
    for(auto iter = Seg.begin(); iter != std::prev(Seg.end()); ++iter)
    {
        double N_i = iter->first;
        double N_f = std::next(iter)->first;
        int count = 0;
        
        while(count == 0)
        {
            double N_m = 0.5 * (N_i + N_f);
            double w_2_i = w_2(N_i, k);
            double w_2_f = w_2(N_f, k);
            double w_2_m = w_2(N_m, k);
            
            //Lin or Exp
            if(w_2_i > 0 and w_2_f > 0)
            {
                //Try lin, log for one segment
                double b_lin_0 = (w_2_f - w_2_i) / (N_f - N_i);
                double a_lin_0 = w_2_i - b_lin_0 * N_i;
                Eigen::Matrix2d Mat_lin_0 = Airy_Mat(a_lin_0, b_lin_0, N_i, N_f);
                double b_log_0 = 0.5 * (log(w_2_f) - log(w_2_i)) / (N_f - N_i);
                double a_log_0 = 0.5 * log(w_2_i) - b_log_0 * N_i;
                Eigen::Matrix2d Mat_log_0 = Bessel_Mat(a_log_0, b_log_0, N_i, N_f);
                
                //Try lin, log for two segment
                double b_lin_1 = (w_2_m - w_2_i) / (N_m - N_i);
                double a_lin_1 = w_2_i - b_lin_1 * N_i;
                double b_lin_2 = (w_2_f - w_2_m) / (N_f - N_m);
                double a_lin_2 = w_2_m - b_lin_2 * N_m;
                Eigen::Matrix2d Mat_lin_1 = Airy_Mat(a_lin_1, b_lin_1, N_i, N_m);
                Eigen::Matrix2d Mat_lin_2 = Airy_Mat(a_lin_2, b_lin_2, N_m, N_f);
                
                double b_log_1 = 0.5 * (log(w_2_m) - log(w_2_i)) / (N_m - N_i);
                double a_log_1 = 0.5 * log(w_2_i) - b_log_1 * N_i;
                double b_log_2 = 0.5 * (log(w_2_f) - log(w_2_m)) / (N_f - N_m);
                double a_log_2 = 0.5 * log(w_2_m) - b_log_2 * N_m;
                Eigen::Matrix2d Mat_log_1 = Bessel_Mat(a_log_1, b_log_1, N_i, N_m);
                Eigen::Matrix2d Mat_log_2 = Bessel_Mat(a_log_2, b_log_2, N_m, N_f);
                
                double Q_lin_1 = abs((Mat_lin_0 * Q_i)[0]);
                double Q_lin_2 = abs((Mat_lin_2 * Mat_lin_1 * Q_i)[0]);
                double err_lin = abs(1 - pow(Q_lin_1 / Q_lin_2, 2));
                
                double Q_log_1 = abs((Mat_log_0 * Q_i)[0]);
                double Q_log_2 = abs((Mat_log_2 * Mat_log_1 * Q_i)[0]);
                double err_log = abs(1 - pow(Q_log_1 / Q_log_2, 2));
                
                if(err_lin < PPS_error or err_log < PPS_error)
                {
                    if(err_lin < err_log)
                        Seg.at(N_i) = Vars(a_lin_0, b_lin_0, 1, Mat_lin_0);
                    else
                        Seg.at(N_i) = Vars(a_log_0, b_log_0, 2, Mat_log_0);
                    
                    Seg.insert(std::pair<double,Vars> (N_f, Null_var));
                    count += 1;
                }
                else
                {
                    N_f = N_m;
                }
            }
            //Lin or -Exp
            else if(w_2_i < 0 and w_2_f < 0)
            {
                //Try lin, log for one segment
                double b_lin_0 = (w_2_f - w_2_i) / (N_f - N_i);
                double a_lin_0 = w_2_i - b_lin_0 * N_i;
                Eigen::Matrix2d Mat_lin_0 = Airy_Mat(a_lin_0, b_lin_0, N_i, N_f);
                double b_log_0 = 0.5 * (log(-w_2_f) - log(-w_2_i)) / (N_f - N_i);
                double a_log_0 = 0.5 * log(-w_2_i) - b_log_0 * N_i;
                Eigen::Matrix2d Mat_log_0 = Modified_Bessel_Mat(a_log_0, b_log_0, N_i, N_f);
                
                //Try lin, log for two segment
                double b_lin_1 = (w_2_m - w_2_i) / (N_m - N_i);
                double a_lin_1 = w_2_i - b_lin_1 * N_i;
                double b_lin_2 = (w_2_f - w_2_m) / (N_f - N_m);
                double a_lin_2 = w_2_m - b_lin_2 * N_m;
                Eigen::Matrix2d Mat_lin_1 = Airy_Mat(a_lin_1, b_lin_1, N_i, N_m);
                Eigen::Matrix2d Mat_lin_2 = Airy_Mat(a_lin_2, b_lin_2, N_m, N_f);
                
                double b_log_1 = 0.5 * (log(-w_2_m) - log(-w_2_i)) / (N_m - N_i);
                double a_log_1 = 0.5 * log(-w_2_i) - b_log_1 * N_i;
                double b_log_2 = 0.5 * (log(-w_2_f) - log(-w_2_m)) / (N_f - N_m);
                double a_log_2 = 0.5 * log(-w_2_m) - b_log_2 * N_m;
                Eigen::Matrix2d Mat_log_1 = Modified_Bessel_Mat(a_log_1, b_log_1, N_i, N_m);
                Eigen::Matrix2d Mat_log_2 = Modified_Bessel_Mat(a_log_2, b_log_2, N_m, N_f);
                
                double Q_lin_1 = abs((Mat_lin_0 * Q_i)[0]);
                double Q_lin_2 = abs((Mat_lin_2 * Mat_lin_1 * Q_i)[0]);
                double err_lin = abs(1 - pow(Q_lin_1 / Q_lin_2, 2));
                
                double Q_log_1 = abs((Mat_log_0 * Q_i)[0]);
                double Q_log_2 = abs((Mat_log_2 * Mat_log_1 * Q_i)[0]);
                double err_log = abs(1 - pow(Q_log_1 / Q_log_2, 2));
                
                if(err_lin < PPS_error or err_log < PPS_error)
                {
                    if(err_lin < err_log)
                        Seg.at(N_i) = Vars(a_lin_0, b_lin_0, 1, Mat_lin_0);
                    else
                        Seg.at(N_i) = Vars(a_log_0, b_log_0, 3, Mat_log_0);
                    
                    Seg.insert(std::pair<double,Vars> (N_f, Null_var));
                    count += 1;
                }
                else
                {
                    N_f = N_m;
                }
            }
            //Lin
            else
            {
                //Try lin one segment
                double b_lin_0 = (w_2_f - w_2_i) / (N_f - N_i);
                double a_lin_0 = w_2_i - b_lin_0 * N_i;
                Eigen::Matrix2d Mat_lin_0 = Airy_Mat(a_lin_0, b_lin_0, N_i, N_f);
                
                //Try lin for two segment
                double b_lin_1 = (w_2_m - w_2_i) / (N_m - N_i);
                double a_lin_1 = w_2_i - b_lin_1 * N_i;
                double b_lin_2 = (w_2_f - w_2_m) / (N_f - N_m);
                double a_lin_2 = w_2_m - b_lin_2 * N_m;
                Eigen::Matrix2d Mat_lin_1 = Airy_Mat(a_lin_1, b_lin_1, N_i, N_m);
                Eigen::Matrix2d Mat_lin_2 = Airy_Mat(a_lin_2, b_lin_2, N_m, N_f);
                
                double Q_lin_1 = abs((Mat_lin_0 * Q_i)[0]);
                double Q_lin_2 = abs((Mat_lin_2 * Mat_lin_1 * Q_i)[0]);
                double err_lin = abs(1 - pow(Q_lin_1 / Q_lin_2, 2));
                
                if(err_lin < PPS_error)
                {
                    Seg.at(N_i) = Vars(a_lin_0, b_lin_0, 1, Mat_lin_0);
                    Seg.insert(std::pair<double,Vars> (N_f, Null_var));
                    count += 1;
                }
                else
                {
                    N_f = N_m;
                }
            }
        }
    }
    
    Eigen::Matrix2d Evolve = Eigen::Matrix2d::Identity();
    for(auto iter = Seg.begin(); iter != std::prev(Seg.end()); ++iter)
        Evolve = iter->second.Mat * Evolve;
    
    Eigen::Vector2cd Q_f = Evolve * Q_i;
    
    return Q_f;
}

double ModeSolver::w_2(double N, double k)
{
    return Bsol.omega_2(N) + k * k * exp(-2 * Bsol.log_aH(N));
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
        
        auto I1 = Bessel_I(0, x1 / p);
        auto K1 = Bessel_K(0, x1 / p);
        auto Ip1 = x1 * Bessel_I(1, x1 / p);
        auto Kp1 = -x1 * Bessel_K(1, x1 / p);
        
        MB0 << I0,  K0,
        Ip0, Kp0;
        
        MB1 << I1,  K1,
        Ip1, Kp1;
        
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
