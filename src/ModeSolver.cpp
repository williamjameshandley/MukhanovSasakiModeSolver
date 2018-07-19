#include "ModeSolver.hpp"

void ModeSolver::Find_Mat()
{
    size_t Nstep = eta_step.size();
    
    std::vector<std::vector<Eigen::MatrixXd>> Mat_vec(Nstep - 1, std::vector<Eigen::MatrixXd>(k.size()));
    Eigen::MatrixXd M(2,2), M_prev;
    
    for(size_t i = 0; i < k.size(); i++)
    {
        M_prev = Eigen::MatrixXd::Identity(2,2);
        
        for(size_t j = Nstep - 2; j != static_cast<size_t>(-1); j--)
        {
            double Ai, Bi, Aip, Bip, Ai0, Bi0, Aip0, Bip0, det;
            
            Ai0 = Airy_Ai((a[j] + b[j] * eta_step[j] - k[i] * k[i]) / std::pow(b[j], 2.0/3.0));
            Bi0 = Airy_Bi((a[j] + b[j] * eta_step[j] - k[i] * k[i]) / std::pow(b[j], 2.0/3.0));
            Aip0 = std::pow(b[j], 1.0/3) * Airy_Aip((a[j] + b[j] * eta_step[j] - k[i] * k[i]) / std::pow(b[j], 2.0/3.0));
            Bip0 = std::pow(b[j], 1.0/3) * Airy_Bip((a[j] + b[j] * eta_step[j] - k[i] * k[i]) / std::pow(b[j], 2.0/3.0));
            
            Ai = Airy_Ai((a[j] + b[j] * eta_step[j+1] - k[i] * k[i]) / std::pow(b[j], 2.0/3.0));
            Bi = Airy_Bi((a[j] + b[j] * eta_step[j+1] - k[i] * k[i]) / std::pow(b[j], 2.0/3.0));
            Aip = std::pow(b[j], 1.0/3) * Airy_Aip((a[j] + b[j] * eta_step[j+1] - k[i] * k[i]) / std::pow(b[j], 2.0/3.0));
            Bip = std::pow(b[j], 1.0/3) * Airy_Bip((a[j] + b[j] * eta_step[j+1] - k[i] * k[i]) / std::pow(b[j], 2.0/3.0));
            
            det = Ai0 * Bip0 - Bi0 * Aip0;
            
            M(0,0) = (Ai*Bip0 - Bi*Aip0) / det;
            M(0,1) = (Bi*Ai0 - Ai*Bi0) / det;
            M(1,0) = (Aip*Bip0 - Bip*Aip0) / det;
            M(1,1) = (Bip*Ai0 - Aip*Bi0) / det;
            
            M_prev = Mat_vec[j][i] = M_prev * M;
        }
        if(i%100 == 0)
        {
            std::cout<<"k = "<<k[i]<<std::endl;
        }
    }
    
    Mat = Mat_vec;
    
}

void ModeSolver::Initial_Conditions(std::string Vacuum, double eta_r)
{
    size_t initial_index = static_cast<size_t>(std::lower_bound(eta_step.begin(), eta_step.end(), eta_r) - eta_step.begin());
    
    std::vector<Eigen::MatrixXd> M_index = Mat[initial_index];
    
    std::vector<std::complex<double>> kin(k.size()), dkin(k.size()), t1(k.size()), t2(k.size());
    
    for(size_t n = 0; n < k.size(); n++)
    {
        kin[n] = 1.0 / sqrt(2 * k[n]);
    }
    
    if(Vacuum == "BD")
    {
        for(size_t n = 0; n < k.size(); n++)
        {
            dkin[n] = -I * k[n] * kin[n];
        }
    }
    else
    {
        std::cout<<"I do not know this initial condition "<<Vacuum<<std::endl;
        exit(0);
    }
    
    if(initial_index == 0)
    {
        for(size_t n = 0; n < k.size(); n++)
        {
            std::complex<double> H10, H20, dH10, dH20, H1, H2, dH1, dH2, det;
            
            H10 = sqrt(eta_r) * Hankel1(0, k[n] * eta_r);
            H20 = sqrt(eta_r) * Hankel2(0, k[n] * eta_r);
            dH10 = 0.5 * Hankel1(0, k[n] * eta_r) / sqrt(eta_r) - k[n] * sqrt(eta_r) * Hankel1(1, k[n] * eta_r);
            dH20 = 0.5 * Hankel2(0, k[n] * eta_r) / sqrt(eta_r) - k[n] * sqrt(eta_r) * Hankel2(1, k[n] * eta_r);
            
            H1 = sqrt(eta_step[0]) * Hankel1(0, k[n] * eta_step[0]);
            H2 = sqrt(eta_step[0]) * Hankel2(0, k[n] * eta_step[0]);
            dH1 = 0.5 * Hankel1(0, k[n] * eta_step[0]) / sqrt(eta_step[0]) - k[n] * sqrt(eta_step[0]) * Hankel1(1, k[n] * eta_step[0]);
            dH2 = 0.5 * Hankel2(0, k[n] * eta_step[0]) / sqrt(eta_step[0]) - k[n] * sqrt(eta_step[0]) * Hankel2(1, k[n] * eta_step[0]);
            
            det = (H10 * dH20 - H20 * dH10);
      
            t1[n] = ((H1 * dH20 - H2 * dH10) * kin[n] + (H2 * H10 - H1 * H20) * dkin[n]) / det;
            t2[n] = ((dH1 * dH20 - dH2 * dH10) * kin[n] + (dH2 * H10 - dH1 * H20) * dkin[n]) / det;
        }
    }
    else if(initial_index != 0)
    {
        for(size_t n = 0; n < k.size(); n++)
        {
            double Ai, Bi, Aip, Bip, Ai0, Bi0, Aip0, Bip0, det;
            
            Ai0 = Airy_Ai((a[initial_index-1] - k[n] * k[n] + b[initial_index-1] * eta_r) / std::pow(b[initial_index-1], 2.0/3));
            Bi0 = Airy_Bi((a[initial_index-1] - k[n] * k[n] + b[initial_index-1] * eta_r) / std::pow(b[initial_index-1], 2.0/3));
            Aip0 = std::pow(b[initial_index-1], 1.0/3) * Airy_Aip((a[initial_index-1] - k[n] * k[n] + b[initial_index-1]*eta_r) / std::pow(b[initial_index-1], 2.0/3));
            Bip0 = std::pow(b[initial_index-1], 1.0/3) * Airy_Bip((a[initial_index-1] - k[n] * k[n] + b[initial_index-1]*eta_r) / std::pow(b[initial_index-1], 2.0/3));
            
            Ai = Airy_Ai((a[initial_index-1] - k[n] * k[n] + b[initial_index-1] * eta_step[initial_index]) / std::pow(b[initial_index-1], 2.0/3));
            Bi = Airy_Bi((a[initial_index-1] - k[n] * k[n] + b[initial_index-1] * eta_step[initial_index]) / std::pow(b[initial_index-1], 2.0/3));
            Aip = std::pow(b[initial_index-1], 1.0/3) * Airy_Aip((a[initial_index-1] - k[n] * k[n] + b[initial_index-1] * eta_step[initial_index]) / std::pow(b[initial_index-1], 2.0/3));
            Bip = std::pow(b[initial_index-1], 1.0/3) * Airy_Bip((a[initial_index-1] - k[n] * k[n] + b[initial_index-1] * eta_step[initial_index]) / std::pow(b[initial_index-1], 2.0/3));
            
            det = Ai0 * Bip0 - Bi0 * Aip0;
            
            t1[n] = ((Ai*Bip0 - Bi*Aip0) * kin[n] + (Bi*Ai0 - Ai*Bi0) * dkin[n]) / det;
            t2[n] = ((Aip*Bip0 - Bip*Aip0) * kin[n] + (Bip*Ai0 - Aip*Bi0) * dkin[n]) / det;
        }
    }
    
    for(size_t n = 0; n < k.size(); n++)
    {
        t1[n] = M_index[n](0,0) * t1[n] + M_index[n](0,1) * t2[n];
        t2[n] = M_index[n](1,0) * t1[n] + M_index[n](1,1) * t2[n];
    }
    
    //MdS
    double v = 1.5;// * sqrt(1 + 8.0 * delta / 9);
    double x = eta_step.back() - eta_end;
    
    for(size_t n = 0; n < k.size(); n++)
    {
        std::complex<double> H1, H2, dH1, dH2;
        
        H1 = I * sqrt(-x) * Hankel1(v, k[n] * x);
        H2 = I * sqrt(-x) * Hankel2(v, k[n] * x);
        dH1 = 0.5 * I * (sqrt(-x) * k[n] * (Hankel1(v-1, k[n] * x) - Hankel1(v+1, k[n] * x)) + Hankel1(v, k[n] * x) / sqrt(-x));
        dH2 = 0.5 * I * (sqrt(-x) * k[n] * (Hankel2(v-1, k[n] * x) - Hankel2(v+1, k[n] * x)) + Hankel2(v, k[n] * x) / sqrt(-x));
        
        c.push_back((dH2 * t1[n] - H2 * t2[n]) / (H1 * dH2 - H2 * dH1));
        d.push_back((-dH1 * t1[n] + H1 * t2[n]) / (H1 * dH2 - H2 * dH1));
        
    }
    
}

std::vector<double> ModeSolver::PPS()
{
    std::complex<double> a11, a12, a21, a22, A_z;
    std::vector<double> PPS(k.size());
    double v = 1.5 * sqrt(1 + 8.0 * delta / 9);
    
    LinearInterpolator<double, double> DZ, Z;
    for(size_t o = 0; o < dz.size(); o++)
    {
        DZ.insert(eta_sol[o], dz[o]);
        Z.insert(eta_sol[o], z[o]);
    }
    
    std::complex<double> x(eta_step.back() - eta_end, 0);
    
    a11 = std::pow(x, (0.5 - v));
    a12 = std::pow(x, (0.5 + v));
    a21 = (0.5 - v) * std::pow(x, (-0.5 - v));
    a22 = (0.5 + v) * std::pow(x, (-0.5 + v));
    
    A_z = (a22 * Z(eta_step.back()) - a12 * DZ(eta_step.back())) / (a11*a22 - a12*a21);
    
    for(size_t n = 0; n < k.size(); n++)
    {
        PPS[n] = (std::pow(k[n], (3.0 - 2.0*v)) / (2.0 * M_PI * M_PI)) * std::pow(std::abs((std::pow(2, v) * Gamma(v) / (A_z * M_PI)) * (c[n] - d[n])), 2);
    }
    
    return PPS;
}
