#include "ModeSolver.hpp"
#include <fstream>

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

Eigen::MatrixXd ModeSolver::lin_step(double w_2_i, double w_2_f, double N_i, double N_f)
{
    double b_lin = (w_2_f - w_2_i) / (N_f - N_i);
    double a_lin = w_2_i - b_lin * N_i;
    return Airy_Mat(a_lin, b_lin, N_i, N_f);
}

Eigen::MatrixXd ModeSolver::pos_exp_step(double w_2_i, double w_2_f, double N_i, double N_f)
{
    double b_log = 0.5 * (log(w_2_f) - log(w_2_i)) / (N_f - N_i);
    double a_log = 0.5 * log(w_2_i) - b_log * N_i;
    return Bessel_Mat(a_log, b_log, N_i, N_f);
}

Eigen::MatrixXd ModeSolver::neg_exp_step(double w_2_i, double w_2_f, double N_i, double N_f)
{
    double b_log = 0.5 * (log(-w_2_f) - log(-w_2_i)) / (N_f - N_i);
    double a_log = 0.5 * log(-w_2_i) - b_log * N_i;
    return Modified_Bessel_Mat(a_log, b_log, N_i, N_f);
}

Eigen::Vector2cd ModeSolver::Evolve(Eigen::Vector2cd Q_i, double k, double N_initial, double& N_final)
{
    std::map<double, Eigen::VectorXcd> Seg;
    
    //Initialize Extrema and end points
    Seg[N_initial] =  Q_i;
    Seg[N_final] =  {};

    for (auto N : Bsol.N_extrema)
        if (N>N_initial and N<N_final)
            Seg[N] = {};

    for (auto i=1;i<99;i++)
        Seg[N_initial + i*(N_final-N_initial)/100] = {};

    
    for(auto iter = Seg.begin(); iter != std::prev(Seg.end()); ++iter)
    {
        double N_i = iter->first;
        double N_f = std::next(iter)->first;
        
        //std::cout << std::log(k) << " " << Bsol.log_aH(N_i) << std::endl;
        if (std::log(k) < std::log(1e-2) +  Bsol.log_aH(N_i))
        {
            N_final = N_i;
            break;
        }
        
        while(true)
        {
            double N_m1 = (2*N_i/3 + N_f/3);
            double N_m2 = (N_i/3 + 2*N_f/3);
            double w_2_i = w_2(N_i, k);
            double w_2_f = w_2(N_f, k);
            double w_2_m1 = w_2(N_m1, k);
            double w_2_m2 = w_2(N_m2, k);

            Eigen::Vector2cd Q_lin_1 = lin_step(w_2_i, w_2_f, N_i, N_f) * iter->second;
            Eigen::VectorXcd Q_lin_2 = lin_step(w_2_m2, w_2_f, N_m2, N_f) * lin_step(w_2_m1, w_2_m2, N_m1, N_m2) * lin_step(w_2_i, w_2_m1, N_i, N_m1) * iter->second;
            double err_lin = std::abs(std::log(abs(Q_lin_2[0]))- std::log(abs(Q_lin_1[0])));
            
            if(w_2_i > 0 and w_2_f > 0)
            {
                Eigen::VectorXcd Q_pos_1 = pos_exp_step(w_2_i, w_2_f, N_i, N_f) * iter->second;
                Eigen::VectorXcd Q_pos_2 = pos_exp_step(w_2_m2, w_2_f, N_m2, N_f) * pos_exp_step(w_2_m1, w_2_m2, N_m1, N_m2) * pos_exp_step(w_2_i, w_2_m1, N_i, N_m1) * iter->second;
                double err_pos = std::abs(std::log(abs(Q_pos_2[0]))- std::log(abs(Q_pos_1[0])));
                
                if(err_lin < PPS_error or err_pos < PPS_error)
                {
                    if (err_lin < err_pos) Seg[N_f] = Q_lin_2; 
                    else                   Seg[N_f] = Q_pos_2; 
                    break;
                }
            }
            else if(w_2_i < 0 and w_2_f < 0)
            {
                Eigen::VectorXcd Q_neg_1 = neg_exp_step(w_2_i, w_2_f, N_i, N_f) * iter->second;
                Eigen::VectorXcd Q_neg_2 = neg_exp_step(w_2_m2, w_2_f, N_m2, N_f) * neg_exp_step(w_2_m1, w_2_m2, N_m1, N_m2) * neg_exp_step(w_2_i, w_2_m1, N_i, N_m1) * iter->second;
                double err_neg = std::abs(std::log(abs(Q_neg_2[0]))- std::log(abs(Q_neg_1[0])));
                
                if(err_lin < PPS_error or err_neg < PPS_error)
                {
                    if (err_lin < err_neg) Seg[N_f] = Q_lin_2;
                    else                   Seg[N_f] = Q_neg_2; 
                    break;
                }
            }
            else
            {
                if(err_lin < PPS_error)
                {
                    Seg[N_f] = Q_lin_2;
                    break;
                }
            }
            Seg[N_m1] = {};
            N_f = N_m1;
        }
    }
    return Seg[N_final];
}

double ModeSolver::w_2(double N, double k)
{
    return Bsol.omega_2(N) + k * k * exp(-2 * Bsol.log_aH(N));
}

Eigen::Matrix2d ModeSolver::Airy_Mat(double a, double b, double N0, double N1)
{
    double p = pow(abs(b), 2.0/3.0);
    double x0 = -((a + b * N0) /p);
    double x1 = -((a + b * N1) /p);

    Eigen::Matrix2d A1, A;
    
    if(x0 < 25 and x1 < 25)
    {
        double Ai0, Bi0, Aip0, Bip0, Ai1, Bi1, Aip1, Bip1;
        Airy(x0, Ai0, Aip0, Bi0, Bip0);
        Airy(x1, Ai1, Aip1, Bi1, Bip1);

        A << (Ai1*Bip0-Aip0*Bi1),          (Ai1*Bi0-Ai0*Bi1)*p/b,
             (Aip0*Bip1-Aip1*Bip0) * b/p,  (Ai0*Bip1-Aip1*Bi0);
        A*=M_PI;
    }
    else
    {
        auto b0 = pow(x0, 0.25), b1 = pow(x1, 0.25);
        auto x = 2.0/3 * (pow(x1, 1.5) - pow(x0, 1.5));
        
        A << std::cosh(x) * b0 / b1,          std::sinh(x) / (-b/p * b0 * b1),
             std::sinh(x) * (-b/p * b0 * b1), std::cosh(x) * b1 / b0;
    }
    
    return A;
}

Eigen::Matrix2d ModeSolver::Bessel_Mat(double a, double b, double N0, double N1)
{
    Eigen::Matrix2d A;
    double p = std::abs(b);
    double x0 = exp(a + b * N0)/std::abs(b);
    double x1 = exp(a + b * N1)/std::abs(b);

    double J0_0 = j0(x0), J1_0 = j1(x0), Y0_0 = y0(x0), Y1_0 = y1(x0);
    double J0_1 = j0(x1), J1_1 = j1(x1), Y0_1 = y0(x1), Y1_1 = y1(x1);

    A <<  (J1_0*Y0_1-J0_1*Y1_0)*x0,        (J0_0*Y0_1-J0_1*Y0_0) / b,
          (J1_1*Y1_0-J1_0*Y1_1)*x0*x1 * b ,(J1_1*Y0_0-J0_0*Y1_1)*x1;
    A *= M_PI/2;
    return A;
}

Eigen::Matrix2d ModeSolver::Modified_Bessel_Mat(double a, double b, double N0, double N1)
{
    double p = b;
    double x0 = exp(a + b * N0);
    double x1 = exp(a + b * N1);
    
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
