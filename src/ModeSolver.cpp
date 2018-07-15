#include "ModeSolver.hpp"

void ModeSolver::Find_Mat()
{
    size_t Nstep = eta_step.size();
    std::vector<std::vector<Eigen::MatrixXd>>
    Mat_vec(k.size(), std::vector<Eigen::MatrixXd>(Nstep - 1));
    
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
        
        Mat_vec[i][eta_step.size() - 2] = M;
        
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
            
            Mat_vec[i][j] = Mat_vec[i][j+1] * M;
        }
        
        std::cout<<"k = "<<k[i]<<std::endl;
    }
    
}
