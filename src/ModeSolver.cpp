#include "ModeSolver.hpp"

ModeSolver::ModeSolver(std::vector<double> kk, 
        std::vector<double> ee, 
        std::vector<double> aa, 
        std::vector<double> bb, 
        double dd, double endend, 
        std::vector<double> zz, 
        std::vector<double> ddzz, 
        std::vector<double> ddddzz, 
        std::vector<double> esol):
            k(kk), eta_step(ee), a(aa), b(bb), delta(dd), eta_end(endend), z(zz), dz(ddzz), ddz(ddddzz), eta_sol(esol), 
            Mat(eta_step.size() - 1, std::vector<Eigen::Matrix2d>(k.size())), c{}, d{}
{}


void ModeSolver::Find_Mat()
{
    for(size_t i = 0; i < k.size(); i++)
    {
        Eigen::Matrix2d M_prev = Eigen::Matrix2d::Identity();
        
        for(size_t j = Mat.size()-1; j != static_cast<size_t>(-1); j--)
        {
            auto Ai0 = Airy_Ai((a[j] + b[j] * eta_step[j] - k[i] * k[i]) / pow(b[j], 2.0/3.0));
            auto Bi0 = Airy_Bi((a[j] + b[j] * eta_step[j] - k[i] * k[i]) / pow(b[j], 2.0/3.0));
            auto Aip0 = pow(b[j], 1.0/3.0) * Airy_Aip((a[j] + b[j] * eta_step[j] - k[i] * k[i]) / pow(b[j], 2.0/3.0));
            auto Bip0 = pow(b[j], 1.0/3.0) * Airy_Bip((a[j] + b[j] * eta_step[j] - k[i] * k[i]) / pow(b[j], 2.0/3.0));
            
            auto Ai = Airy_Ai((a[j] + b[j] * eta_step[j+1] - k[i] * k[i]) / pow(b[j], 2.0/3.0));
            auto Bi = Airy_Bi((a[j] + b[j] * eta_step[j+1] - k[i] * k[i]) / pow(b[j], 2.0/3.0));
            auto Aip = pow(b[j], 1.0/3.0) * Airy_Aip((a[j] + b[j] * eta_step[j+1] - k[i] * k[i]) / pow(b[j], 2.0/3.0));
            auto Bip = pow(b[j], 1.0/3.0) * Airy_Bip((a[j] + b[j] * eta_step[j+1] - k[i] * k[i]) / pow(b[j], 2.0/3.0));
            
            auto det = Ai0 * Bip0 - Bi0 * Aip0;
            
            Eigen::Matrix2d M;
            M << (Ai*Bip0 - Bi*Aip0) / det,   (Bi*Ai0 - Ai*Bi0) / det,
                 (Aip*Bip0 - Bip*Aip0) / det, (Bip*Ai0 - Aip*Bi0) / det; 
            
            M_prev = Mat[j][i] = M_prev * M;
        }
        if(i%100 == 0)
            std::cout<<"k = "<<k[i]<<std::endl;
    }
    
}

void ModeSolver::Initial_Conditions(std::string Vacuum, double eta_r)
{
    
    std::vector<std::complex<double>> kin(k.size()), dkin(k.size()), t1(k.size()), t2(k.size());
    for(size_t n = 0; n < k.size(); n++)
        kin[n] = 1.0 / sqrt(2 * k[n]);
    
    if(Vacuum == "BD")
        for(size_t n = 0; n < k.size(); n++)
            dkin[n] = -I * k[n] * kin[n];
    else
        throw std::runtime_error("Initial conditions unknown");
    
    auto initial_index = static_cast<size_t>(std::lower_bound(eta_step.begin(), eta_step.end(), eta_r) - eta_step.begin());

    if(initial_index == 0)
        for(size_t n = 0; n < k.size(); n++)
        {
            auto H10 = sqrt(eta_r) * Hankel1(0, k[n] * eta_r);
            auto H20 = sqrt(eta_r) * Hankel2(0, k[n] * eta_r);
            auto dH10 = 0.5 * Hankel1(0, k[n] * eta_r) / sqrt(eta_r) - k[n] * sqrt(eta_r) * Hankel1(1, k[n] * eta_r);
            auto dH20 = 0.5 * Hankel2(0, k[n] * eta_r) / sqrt(eta_r) - k[n] * sqrt(eta_r) * Hankel2(1, k[n] * eta_r);
            
            auto H1 = sqrt(eta_step[0]) * Hankel1(0, k[n] * eta_step[0]);
            auto H2 = sqrt(eta_step[0]) * Hankel2(0, k[n] * eta_step[0]);
            auto dH1 = 0.5 * Hankel1(0, k[n] * eta_step[0]) / sqrt(eta_step[0]) - k[n] * sqrt(eta_step[0]) * Hankel1(1, k[n] * eta_step[0]);
            auto dH2 = 0.5 * Hankel2(0, k[n] * eta_step[0]) / sqrt(eta_step[0]) - k[n] * sqrt(eta_step[0]) * Hankel2(1, k[n] * eta_step[0]);
            
            auto det = (H10 * dH20 - H20 * dH10);
      
            t1[n] = ((H1 * dH20 - H2 * dH10) * kin[n] + (H2 * H10 - H1 * H20) * dkin[n]) / det;
            t2[n] = ((dH1 * dH20 - dH2 * dH10) * kin[n] + (dH2 * H10 - dH1 * H20) * dkin[n]) / det;
        }
    else
        for(size_t n = 0; n < k.size(); n++)
        {
            auto Ai0 = Airy_Ai((a[initial_index-1] + b[initial_index-1] * eta_r - k[n] * k[n]) / pow(b[initial_index-1], 2.0/3.0));
            auto Bi0 = Airy_Bi((a[initial_index-1] + b[initial_index-1] * eta_r - k[n] * k[n]) / pow(b[initial_index-1], 2.0/3.0));
            auto Aip0 = pow(b[initial_index-1], 1.0/3.0) * Airy_Aip((a[initial_index-1] + b[initial_index-1] * eta_r - k[n] * k[n]) / pow(b[initial_index-1], 2.0/3.0));
            auto Bip0 = pow(b[initial_index-1], 1.0/3.0) * Airy_Bip((a[initial_index-1] + b[initial_index-1] * eta_r - k[n] * k[n]) / pow(b[initial_index-1], 2.0/3.0));
            
            auto Ai = Airy_Ai((a[initial_index-1] + b[initial_index-1] * eta_step[initial_index] - k[n] * k[n]) / pow(b[initial_index-1], 2.0/3));
            auto Bi = Airy_Bi((a[initial_index-1] + b[initial_index-1] * eta_step[initial_index] - k[n] * k[n]) / pow(b[initial_index-1], 2.0/3));
            auto Aip = pow(b[initial_index-1], 1.0/3.0) * Airy_Aip((a[initial_index-1] + b[initial_index-1] * eta_step[initial_index] - k[n] * k[n]) / pow(b[initial_index-1], 2.0/3));
            auto Bip = pow(b[initial_index-1], 1.0/3.0) * Airy_Bip((a[initial_index-1] + b[initial_index-1] * eta_step[initial_index] - k[n] * k[n]) / pow(b[initial_index-1], 2.0/3));
            
            auto det = Ai0 * Bip0 - Bi0 * Aip0;
            
            t1[n] = ((Ai * Bip0 - Bi * Aip0) * kin[n] + (Bi * Ai0 - Ai * Bi0) * dkin[n]) / det;
            t2[n] = ((Aip * Bip0 - Bip * Aip0) * kin[n] + (Bip * Ai0 - Aip * Bi0) * dkin[n]) / det;
        }
    
    
    std::vector<Eigen::Matrix2d> M_index = Mat[initial_index];
    for(size_t n = 0; n < k.size(); n++)
    {
        auto temp1 = t1[n];
        auto temp2 = t2[n];
        t1[n] = M_index[n](0,0) * temp1 + M_index[n](0,1) * temp2;
        t2[n] = M_index[n](1,0) * temp1 + M_index[n](1,1) * temp2;
    }
    
    //MdS
    double v = 1.5 * sqrt(1 + 8.0 * delta / 9);
    double x = eta_step.back() - eta_end;
    
    for(size_t n = 0; n < k.size(); n++)
    {
        auto H1 = I * sqrt(-x) * Hankel1(v, k[n] * x);
        auto H2 = I * sqrt(-x) * Hankel2(v, k[n] * x);
        auto dH1 = 0.5 * I * (sqrt(-x) * k[n] * (Hankel1(v-1, k[n] * x) - Hankel1(v+1, k[n] * x)) - Hankel1(v, k[n] * x) / sqrt(-x));
        auto dH2 = 0.5 * I * (sqrt(-x) * k[n] * (Hankel2(v-1, k[n] * x) - Hankel2(v+1, k[n] * x)) - Hankel2(v, k[n] * x) / sqrt(-x));
        
        c.push_back((dH2 * t1[n] - H2 * t2[n]) / (H1 * dH2 - H2 * dH1));
        d.push_back((-dH1 * t1[n] + H1 * t2[n]) / (H1 * dH2 - H2 * dH1));
        
    }
    
}

std::vector<double> ModeSolver::PPS()
{
    
    LinearInterpolator<double, double> DZ, Z;
    for(size_t o = 0; o < dz.size(); o++)
    {
        DZ.insert(eta_sol[o], dz[o]);
        Z.insert(eta_sol[o], z[o]);
    }
    
    std::complex<double> x(eta_step.back() - eta_end, 0);
    
    auto v = 1.5 * sqrt(1 + 8.0 * delta / 9);
    auto a11 = std::pow(x, (0.5 - v));
    auto a12 = std::pow(x, (0.5 + v));
    auto a21 = (0.5 - v) * std::pow(x, (-0.5 - v));
    auto a22 = (0.5 + v) * std::pow(x, (-0.5 + v));
    
    auto A_z = (a22 * Z(eta_step.back()) - a12 * DZ(eta_step.back())) / (a11*a22 - a12*a21);
    
    std::vector<double> PPS(k.size());
    for(size_t n = 0; n < k.size(); n++)
        PPS[n] = (std::pow(k[n], (3.0 - 2.0*v)) / (2.0 * M_PI * M_PI)) * std::pow(std::abs((std::pow(2, v) * Gamma(v) / (A_z * M_PI)) * (c[n] - d[n])), 2);
    
    return PPS;
}
