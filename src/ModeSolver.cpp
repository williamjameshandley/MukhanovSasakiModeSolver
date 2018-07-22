#include "ModeSolver.hpp"

ModeSolver::ModeSolver(
        std::vector<double> ee, 
        std::vector<double> aa, 
        std::vector<double> bb, 
        double dd, double endend, 
        std::vector<double> zz, 
        std::vector<double> ddzz, 
        std::vector<double> ddddzz, 
        std::vector<double> esol):
            eta_step(ee), a(aa), b(bb), delta(dd), eta_end(endend), z(zz), dz(ddzz), ddz(ddddzz), eta_sol(esol), Mat{}, c{}, d{}
{}


void ModeSolver::Find_Mat(double k)
{
    Mat = Eigen::Matrix2d::Identity();
    
    for(size_t j = eta_step.size() - 2; j != initial_index - 1; j--)
    {
        auto Ai0 = Airy_Ai((a[j] + b[j] * eta_step[j] - k * k) / pow(b[j], 2.0/3.0));
        auto Bi0 = Airy_Bi((a[j] + b[j] * eta_step[j] - k * k) / pow(b[j], 2.0/3.0));
        auto Aip0 = pow(b[j], 1.0/3.0) * Airy_Aip((a[j] + b[j] * eta_step[j] - k * k) / pow(b[j], 2.0/3.0));
        auto Bip0 = pow(b[j], 1.0/3.0) * Airy_Bip((a[j] + b[j] * eta_step[j] - k * k) / pow(b[j], 2.0/3.0));
        
        auto Ai = Airy_Ai((a[j] + b[j] * eta_step[j+1] - k * k) / pow(b[j], 2.0/3.0));
        auto Bi = Airy_Bi((a[j] + b[j] * eta_step[j+1] - k * k) / pow(b[j], 2.0/3.0));
        auto Aip = pow(b[j], 1.0/3.0) * Airy_Aip((a[j] + b[j] * eta_step[j+1] - k * k) / pow(b[j], 2.0/3.0));
        auto Bip = pow(b[j], 1.0/3.0) * Airy_Bip((a[j] + b[j] * eta_step[j+1] - k * k) / pow(b[j], 2.0/3.0));
        
        auto det = Ai0 * Bip0 - Bi0 * Aip0;
        
        Eigen::Matrix2d M;
        M << (Ai*Bip0 - Bi*Aip0) / det,   (Bi*Ai0 - Ai*Bi0) / det,
             (Aip*Bip0 - Bip*Aip0) / det, (Bip*Ai0 - Aip*Bi0) / det;
        
        Mat = Mat * M;
    }
}

void ModeSolver::Initial_Conditions(std::string Vac, double eta_rr)
{
    Vacuum = Vac;
    eta_r = eta_rr;
    initial_index = static_cast<size_t>(std::lower_bound(eta_step.begin(), eta_step.end(), eta_r) - eta_step.begin());
    
    for(size_t o = 0; o < dz.size(); o++)
    {
        DDZ.insert(eta_sol[o], ddz[o]);
        DZ.insert(eta_sol[o], dz[o]);
        Z.insert(eta_sol[o], z[o]);
    }
    
}

void ModeSolver::Match(double k)
{
    std::complex<double> kin, dkin, t1, t2;
    
    kin = 1.0 / sqrt(2 * k);
    
    if(Vacuum == "BD")
        dkin = -I * k * kin;
    else if(Vacuum == "HD")
    {
        dkin = -I * std::pow(k * k - DDZ(eta_r), 0.5) * kin;
    }
    else if(Vacuum == "RSET")
    {
        dkin = (-I * k + DZ(eta_r)/Z(eta_r)) * kin;
    }
    else
        throw std::runtime_error("Initial conditions unknown");

    if(initial_index == 0)
    {
        auto H10 = sqrt(eta_r) * Hankel1(0, k * eta_r);
        auto H20 = sqrt(eta_r) * Hankel2(0, k * eta_r);
        auto dH10 = 0.5 * Hankel1(0, k * eta_r) / sqrt(eta_r) - k * sqrt(eta_r) * Hankel1(1, k * eta_r);
        auto dH20 = 0.5 * Hankel2(0, k * eta_r) / sqrt(eta_r) - k * sqrt(eta_r) * Hankel2(1, k * eta_r);
    
        auto H1 = sqrt(eta_step[0]) * Hankel1(0, k * eta_step[0]);
        auto H2 = sqrt(eta_step[0]) * Hankel2(0, k * eta_step[0]);
        auto dH1 = 0.5 * Hankel1(0, k * eta_step[0]) / sqrt(eta_step[0]) - k * sqrt(eta_step[0]) * Hankel1(1, k * eta_step[0]);
        auto dH2 = 0.5 * Hankel2(0, k * eta_step[0]) / sqrt(eta_step[0]) - k * sqrt(eta_step[0]) * Hankel2(1, k * eta_step[0]);
    
        auto det = (H10 * dH20 - H20 * dH10);
  
        t1 = ((H1 * dH20 - H2 * dH10) * kin + (H2 * H10 - H1 * H20) * dkin) / det;
        t2 = ((dH1 * dH20 - dH2 * dH10) * kin + (dH2 * H10 - dH1 * H20) * dkin) / det;
    }
    else
    {
        auto Ai0 = Airy_Ai((a[initial_index-1] + b[initial_index-1] * eta_r - k * k) / pow(b[initial_index-1], 2.0/3.0));
        auto Bi0 = Airy_Bi((a[initial_index-1] + b[initial_index-1] * eta_r - k * k) / pow(b[initial_index-1], 2.0/3.0));
        auto Aip0 = pow(b[initial_index-1], 1.0/3.0) * Airy_Aip((a[initial_index-1] + b[initial_index-1] * eta_r - k * k) / pow(b[initial_index-1], 2.0/3.0));
        auto Bip0 = pow(b[initial_index-1], 1.0/3.0) * Airy_Bip((a[initial_index-1] + b[initial_index-1] * eta_r - k * k) / pow(b[initial_index-1], 2.0/3.0));
    
        auto Ai = Airy_Ai((a[initial_index-1] + b[initial_index-1] * eta_step[initial_index] - k * k) / pow(b[initial_index-1], 2.0/3));
        auto Bi = Airy_Bi((a[initial_index-1] + b[initial_index-1] * eta_step[initial_index] - k * k) / pow(b[initial_index-1], 2.0/3));
        auto Aip = pow(b[initial_index-1], 1.0/3.0) * Airy_Aip((a[initial_index-1] + b[initial_index-1] * eta_step[initial_index] - k * k) / pow(b[initial_index-1], 2.0/3));
        auto Bip = pow(b[initial_index-1], 1.0/3.0) * Airy_Bip((a[initial_index-1] + b[initial_index-1] * eta_step[initial_index] - k * k) / pow(b[initial_index-1], 2.0/3));
    
        auto det = Ai0 * Bip0 - Bi0 * Aip0;
    
        t1 = ((Ai * Bip0 - Bi * Aip0) * kin + (Bi * Ai0 - Ai * Bi0) * dkin) / det;
        t2 = ((Aip * Bip0 - Bip * Aip0) * kin + (Bip * Ai0 - Aip * Bi0) * dkin) / det;
    }
    
    
    Find_Mat(k);
    auto temp1 = t1;
    auto temp2 = t2;
    t1 = Mat(0,0) * temp1 + Mat(0,1) * temp2;
    t2 = Mat(1,0) * temp1 + Mat(1,1) * temp2;
    
    //MdS
    double v = 1.5 * sqrt(1 + 8.0 * delta / 9);
    double x = eta_step.back() - eta_end;

    auto H1 = I * sqrt(-x) * Hankel1(v, k * x);
    auto H2 = I * sqrt(-x) * Hankel2(v, k * x);
    auto dH1 = 0.5 * I * (sqrt(-x) * k * (Hankel1(v-1, k * x) - Hankel1(v+1, k * x)) - Hankel1(v, k * x) / sqrt(-x));
    auto dH2 = 0.5 * I * (sqrt(-x) * k * (Hankel2(v-1, k * x) - Hankel2(v+1, k * x)) - Hankel2(v, k * x) / sqrt(-x));
    
    c = ((dH2 * t1 - H2 * t2) / (H1 * dH2 - H2 * dH1));
    d = ((-dH1 * t1 + H1 * t2) / (H1 * dH2 - H2 * dH1));
    
    
}

double ModeSolver::PPS(double k)
{
    Match(k);
    
    std::complex<double> x(eta_step.back() - eta_end, 0);
    
    auto v = 1.5 * sqrt(1 + 8.0 * delta / 9);
    auto a11 = std::pow(x, (0.5 - v));
    auto a12 = std::pow(x, (0.5 + v));
    auto a21 = (0.5 - v) * std::pow(x, (-0.5 - v));
    auto a22 = (0.5 + v) * std::pow(x, (-0.5 + v));
    
    auto A_z = (a22 * Z(eta_step.back()) - a12 * DZ(eta_step.back())) / (a11*a22 - a12*a21);
    
    return (std::pow(k, (3.0 - 2.0*v)) / (2.0 * M_PI * M_PI)) * std::pow(std::abs((std::pow(2, v) * Gamma(v) / (A_z * M_PI)) * (c - d)), 2);

}
