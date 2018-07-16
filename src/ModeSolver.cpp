#include "ModeSolver.hpp"

void ModeSolver::Find_Mat()
{
    size_t Nstep = eta_step.size();
    std::vector<std::vector<Eigen::MatrixXd>>
    Mat_vec(Nstep - 1, std::vector<Eigen::MatrixXd>(k.size()));
    
    Eigen::MatrixXd M(2,2);
    
    double Ai, Bi, Aip, Bip, Ai0, Bi0, Aip0, Bip0, det;
    
    for(size_t i = 0; i < k.size(); i++)
    {
        Ai0 = boost::math::airy_ai((a[Nstep - 2] - k[i] * k[i] + b[Nstep - 2] * eta_step[Nstep - 2]) / std::pow(b[Nstep - 2], 2.0/3));
        Bi0 = boost::math::airy_bi((a[Nstep - 2] - k[i] * k[i] + b[Nstep - 2] * eta_step[Nstep - 2]) / std::pow(b[Nstep - 2], 2.0/3));
        Aip0 = std::pow(b[Nstep - 2], 1.0/3) * boost::math::airy_ai_prime((a[Nstep - 2] - k[i] * k[i] + b[Nstep - 2] * eta_step[Nstep - 2]) / std::pow(b[Nstep - 2], 2.0/3));
        Bip0 = std::pow(b[Nstep - 2], 1.0/3) * boost::math::airy_bi_prime((a[Nstep - 2] - k[i] * k[i] + b[Nstep - 2] * eta_step[Nstep - 2]) / std::pow(b[Nstep - 2], 2.0/3));
        
        Ai = boost::math::airy_ai((a[Nstep - 2] - k[i] * k[i] + b[Nstep - 2] * eta_step[Nstep - 1]) / std::pow(b[Nstep - 2], 2.0/3));
        Bi = boost::math::airy_bi((a[Nstep - 2] - k[i] * k[i] + b[Nstep - 2] * eta_step[Nstep - 1]) / std::pow(b[Nstep - 2], 2.0/3));
        Aip = std::pow(b[Nstep - 2], 1.0/3) * boost::math::airy_ai_prime((a[Nstep - 2] - k[i] * k[i] + b[Nstep - 2] * eta_step[Nstep - 1]) / std::pow(b[Nstep - 2], 2.0/3));
        Bip = std::pow(b[Nstep - 2], 1.0/3) * boost::math::airy_bi_prime((a[Nstep - 2] - k[i] * k[i] + b[Nstep - 2] * eta_step[Nstep - 1]) / std::pow(b[Nstep - 2], 2.0/3));
        
        det = Ai0 * Bip0 + Bi0 * Aip0;
        
        M(0,0) = (Ai*Bip0 - Bi*Aip0) / det;
        M(0,1) = (Bi*Ai0 - Ai*Bi0) / det;
        M(1,0) = (Aip*Bip0 - Bip*Aip0) / det;
        M(1,1) = (Bip*Ai0 - Aip*Bi0) / det;
        
        Mat_vec[eta_step.size() - 2][i] = M;
        
        for(size_t j = Nstep - 3; j != static_cast<size_t>(-1); j--)
        {
            Ai0 = boost::math::airy_ai((a[j] - k[i] * k[i] + b[j] * eta_step[j]) / std::pow(b[j], 2.0/3));
            Bi0 = boost::math::airy_bi((a[j] - k[i] * k[i] + b[j] * eta_step[j]) / std::pow(b[j], 2.0/3));
            Aip0 = std::pow(b[j], 1.0/3) * boost::math::airy_ai_prime((a[j] - k[i] * k[i] + b[j] * eta_step[j]) / std::pow(b[j], 2.0/3));
            Bip0 = std::pow(b[j], 1.0/3) * boost::math::airy_bi_prime((a[j] - k[i] * k[i] + b[j] * eta_step[j]) / std::pow(b[j], 2.0/3));
            
            Ai = boost::math::airy_ai((a[j] - k[i] * k[i] + b[j] * eta_step[j+1]) / std::pow(b[j], 2.0/3));
            Bi = boost::math::airy_bi((a[j] - k[i] * k[i] + b[j] * eta_step[j+1]) / std::pow(b[j], 2.0/3));
            Aip = std::pow(b[j], 1.0/3) * boost::math::airy_ai_prime((a[j] - k[i] * k[i] + b[j] * eta_step[j+1]) / std::pow(b[j], 2.0/3));
            Bip = std::pow(b[j], 1.0/3) * boost::math::airy_bi_prime((a[j] - k[i] * k[i] + b[j] * eta_step[j+1]) / std::pow(b[j], 2.0/3));
            
            det = Ai0 * Bip0 + Bi0 * Aip0;
            
            M(0,0) = (Ai*Bip0 - Bi*Aip0) / det;
            M(0,1) = (Bi*Ai0 - Ai*Bi0) / det;
            M(1,0) = (Aip*Bip0 - Bip*Aip0) / det;
            M(1,1) = (Bip*Ai0 - Aip*Bi0) / det;
            
            Mat_vec[j][i] = Mat_vec[j+1][i] * M;
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
    size_t i = static_cast<size_t>(std::lower_bound(eta_step.begin(), eta_step.end(), eta_r) - eta_step.begin());
    
    std::vector<Eigen::MatrixXd> M = Mat[i];
    
    std::vector<double> kin;
    std::vector<std::complex<double>> dkin;
    
    std::vector<std::complex<double>> t1(k.size()), t2(k.size()), c(k.size()), d(k.size());
    
    for(size_t n = 0; n < k.size(); n++)
    {
        kin.push_back(1.0 / sqrt(2 * k[n]));
    }
    
    if(Vacuum == "BD")
    {
        for(size_t n = 0; n < k.size(); n++)
        {
            dkin.push_back((0, -k[n] * kin[n]));
        }
    }
    else if(Vacuum == "HD")
    {
        LinearInterpolator<double, double> DDZ;
        for(size_t o = 0; o < ddz.size(); o++)
        {
            DDZ.insert(eta_sol[o], ddz[o]);
        }
        for(size_t n = 0; n < k.size(); n++)
        {
            dkin.push_back((0, -sqrt(k[n] * k[n] - DDZ(eta_r)) * kin[n]));
        }
    }
    else if(Vacuum == "RSET")
    {
        LinearInterpolator<double, double> DZ, Z;
        for(size_t o = 0; o < dz.size(); o++)
        {
            DZ.insert(eta_sol[o], dz[o]);
            Z.insert(eta_sol[o], z[o]);
        }
        for(size_t n = 0; n < k.size(); n++)
        {
            dkin.push_back((DZ(eta_r) * kin[n] / Z(eta_r), -k[n] * kin[n]));
        }
    }
    else
    {
        std::cout<<"I do not know this initial condition "<<Vacuum<<std::endl;
        exit(0);
    }
    
    
    if(i == 0)
    {
        double J0, Y0, dJ0, dY0, J, Y, dJ, dY;
        
        for(size_t n = 0; n < k.size(); n++)
        {
            J0 = sqrt(eta_r) * boost::math::cyl_bessel_j(0, k[n] * eta_r);
            Y0 = sqrt(eta_r) * boost::math::cyl_neumann(0, k[n] * eta_r);
            dJ0 = 0.5 * boost::math::cyl_bessel_j(0, k[n] * eta_r) / sqrt(eta_r) - k[n] * sqrt(eta_r) * boost::math::cyl_bessel_j(1, k[n] * eta_r);
            dY0 = 0.5 * boost::math::cyl_neumann(0, k[n] * eta_r) / sqrt(eta_r) - k[n] * sqrt(eta_r) * boost::math::cyl_neumann(1, k[n] * eta_r);
            
            J = sqrt(eta_step[0]) * boost::math::cyl_bessel_j(0, k[n] * eta_step[0]);
            Y = sqrt(eta_step[0]) * boost::math::cyl_neumann(0, k[n] * eta_step[0]);
            dJ = 0.5 * boost::math::cyl_bessel_j(0, k[n] * eta_step[0]) / sqrt(eta_step[0]) - k[n] * sqrt(eta_step[0]) * boost::math::cyl_bessel_j(1, k[n] * eta_step[0]);
            dY = 0.5 * boost::math::cyl_neumann(0, k[n] * eta_step[0]) / sqrt(eta_step[0]) - k[n] * sqrt(eta_step[0]) * boost::math::cyl_neumann(1, k[n] * eta_step[0]);
            
            t1[n] = ((J*dY0 - Y*dJ0)*kin[i] + (Y*J0 - J*Y0)*dkin[i]) / (J0 * dY0 - Y0 * dJ0);
            t2[n] = ((dJ*dY0 - dY*dJ0)*kin[i] + (dY*J0 - dJ*Y0)*dkin[i]) / (J0 * dY0 - Y0 * dJ0);
        }
    }
    else if(i != 0)
    {
        double Ai, Bi, Aip, Bip, Ai0, Bi0, Aip0, Bip0;
        
        for(size_t n = 0; n < k.size(); n++)
        {
            Ai0 = boost::math::airy_ai((a[i-1] - k[n] * k[n] + b[i-1]*eta_r) / std::pow(b[i-1], 2.0/3));
            Bi0 = boost::math::airy_bi((a[i-1] - k[n] * k[n] + b[i-1]*eta_r) / std::pow(b[i-1], 2.0/3));
            Aip0 = std::pow(b[i-1], 1.0/3) * boost::math::airy_ai_prime((a[i-1] - k[n] * k[n] + b[i-1]*eta_r) / std::pow(b[i-1], 2.0/3));
            Bip0 = std::pow(b[i-1], 1.0/3) * boost::math::airy_bi_prime((a[i-1] - k[n] * k[n] + b[i-1]*eta_r) / std::pow(b[i-1], 2.0/3));
            
            Ai = boost::math::airy_ai((a[i-1] - k[n] * k[n] + b[i-1]*eta_step[i]) / std::pow(b[i-1], 2.0/3));
            Bi = boost::math::airy_bi((a[i-1] - k[n] * k[n] + b[i-1]*eta_step[i]) / std::pow(b[i-1], 2.0/3));
            Aip = std::pow(b[i-1], 1.0/3) * boost::math::airy_ai_prime((a[i-1] - k[n] * k[n] + b[i-1]*eta_step[i]) / std::pow(b[i-1], 2.0/3));
            Bip = std::pow(b[i-1], 1.0/3) * boost::math::airy_bi_prime((a[i-1] - k[n] * k[n] + b[i-1]*eta_step[i]) / std::pow(b[i-1], 2.0/3));
            
            t1[n] = ((Ai*Bip0 - Bi*Aip0)*kin[n] + (Bi*Ai0 - Ai*Bi0)*dkin[n]) / (Ai0*Bip0 - Bi0*Aip0);
            t2[n] = ((Aip*Bip0 - Bip*Aip0)*kin[n] + (Bip*Ai0 - Aip*Bi0)*dkin[n]) / (Ai0*Bip0 - Bi0*Aip0);
        }
    }
    
    for(size_t n = 0; n < k.size(); n++)
    {
        c[n] = M[n](0,0) * t1[n] + M[n](0,1) * t2[n];
        d[n] = M[n](1,0) * t1[n] + M[n](1,1) * t2[n];
        t1[n] = c[n];
        t2[n] = d[n];
    }
    
    
    //MdS
    double v = 1.5 * sqrt(1 + 8.0 * delta / 9);
    std::complex<double> x(eta_step.back() - eta_end,0), BJ, BY, dBJ, dBY;
    
    for(size_t n = 0; n < k.size(); n++)
    {
        BJ = sqrt(x) * boost::math::cyl_bessel_j(v,k[n]*(eta_step.back() - eta_end));
        BY = sqrt(x) * boost::math::cyl_neumann(v,k[n]*(eta_step.back() - eta_end));
        dBJ = 0.5 * (sqrt(x) * k[n] * (boost::math::cyl_bessel_j(v-1,k[n]*(eta_step.back() - eta_end)) - boost::math::cyl_bessel_j(v+1,k[n]*(eta_step.back() - eta_end))) + boost::math::cyl_bessel_j(v,k[n]*(eta_step.back() - eta_end)) / sqrt(x));
        dBY = 0.5 * (sqrt(x) * k[n] * (boost::math::cyl_neumann(v-1,k[n]*(eta_step.back() - eta_end)) - boost::math::cyl_neumann(v+1,k[n]*(eta_step.back() - eta_end))) + boost::math::cyl_neumann(v,k[n]*(eta_step.back() - eta_end)) / sqrt(x));
        
        c[n] = (dBY*t1[n] - BY*t2[n]) / (BJ*dBY - BY*dBJ);
        d[n] = (-dBJ*t1[n] + BJ*t2[n]) / (BJ*dBY - BY*dBJ);
    }
    
    C = c;
    D = d;
    
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
    
    a11 = std::pow((eta_step.back() - eta_end, 0), (0.5 - v));
    a12 = std::pow((eta_step.back() - eta_end, 0), (0.5 + v));
    a21 = (0.5 - v) * std::pow((eta_step.back() - eta_end, 0), (-0.5 - v));
    a22 = (0.5 + v) * std::pow((eta_step.back() - eta_end, 0), (-0.5 + v));
    
    A_z = (a22 * Z(eta_step.back()) - a12 * DZ(eta_step.back())) / (a11*a22 - a12*a21);
    
    for(size_t n = 0; n < k.size(); n++)
    {
        PPS[n] = (std::pow(k[n], (3.-2*v)) / (2 * M_PI * M_PI)) * std::pow(std::abs((std::pow(2, v) * tgamma(v) / (A_z * M_PI)) * (C[n] - D[n])), 2);
    }
    
    return PPS;
}
