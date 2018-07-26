#include "ModeSolver.hpp"

ModeSolver::ModeSolver(BackgroundSolution _Bsol, TransitionsSolution _Tsol):
    Bsol{_Bsol}, Tsol{_Tsol}, eta_r{0.5}, vacuum{BD}, initial_index{0}, DDZ{}, DZ{}, Z{} 
{


}


Eigen::Matrix2d ModeSolver::Mat(double k)
{
    Eigen::Matrix2d Mat = Eigen::Matrix2d::Identity();
    
    for(size_t j = Tsol.eta_step.size() - 2; j != initial_index - 1; j--)
    {
        double p = pow(Tsol.b[j], 1.0/3.0);
        double x = ((Tsol.a[j] + Tsol.b[j] * Tsol.eta_step[j] - k * k) /p/p);
        Eigen::Matrix2d A0 = A(x,p);

        x = ((Tsol.a[j] + Tsol.b[j] * Tsol.eta_step[j+1] - k * k) /p/p);
        Eigen::Matrix2d A1 = A(x,p);

        Mat *=  A1 * A0.inverse();
    }
    return Mat;
}

void ModeSolver::Initial_Conditions(VacuumChoice vac, double eta_rr)
{
    vacuum = vac;
    eta_r = eta_rr;
    initial_index = static_cast<size_t>(std::lower_bound(Tsol.eta_step.begin(), Tsol.eta_step.end(), eta_r) - Tsol.eta_step.begin());
    
    for(size_t o = 0; o < Bsol.dz.size(); o++)
    {
        DDZ.insert(Bsol.eta[o], Bsol.ddz[o]);
        DZ.insert(Bsol.eta[o], Bsol.dz[o]);
        Z.insert(Bsol.eta[o], Bsol.z[o]);
    }
    
}

Eigen::Vector2cd ModeSolver::Match(double k)
{
    Eigen::Vector2cd kin;
    
    kin[0] = 1.0 / sqrt(2 * k);
    
    if(vacuum == BD)
        kin[1] = -I * k * kin[0];
    else if(vacuum == HD)
        kin[1] = -I * std::pow(k * k - DDZ(eta_r), 0.5) * kin[0];
    else if(vacuum == RST)
        kin[1] = (-I * k + DZ(eta_r)/Z(eta_r)) * kin[0];
    else
        throw std::runtime_error("Initial conditions unknown");
    
    double v = 1.5 * sqrt(1 + 8.0 * Tsol.delta / 9);
    double x = Tsol.eta_step.back() - Bsol.eta.back();
    Eigen::Matrix2cd Hi = H(x,k,v);

    if(initial_index == 0)
    {
        Eigen::Matrix2cd H1 = H(eta_r, k); 
        Eigen::Matrix2cd H0 = H(Tsol.eta_step[0], k);

        return Hi.inverse() * Mat(k) * H1 * H0.inverse() * kin;
    }
    else
    {
        double p = pow(Tsol.b[initial_index-1], 1.0/3.0);
        x = ((Tsol.a[initial_index-1] + Tsol.b[initial_index-1] * eta_r - k * k) /p/p);
        Eigen::Matrix2d A0 = A(x,p);

        x = ((Tsol.a[initial_index-1] + Tsol.b[initial_index-1] * Tsol.eta_step[initial_index] - k * k) /p/p);
        Eigen::Matrix2d A1 = A(x,p);

        return Hi.inverse() * Mat(k) * A1 * A0.inverse() * kin;
    }
}

double ModeSolver::PPS(double k)
{
    Eigen::Vector2cd cd = Match(k);
    auto c = cd[0];
    auto d = cd[1];
    
    auto x = Tsol.eta_step.back() - Bsol.eta.back();
    auto v = 1.5 * sqrt(1 + 8.0 * Tsol.delta / 9);
    Eigen::Matrix2cd a_ = a(x, v);
    Eigen::Vector2d z;
    z << Z(Tsol.eta_step.back()), DZ(Tsol.eta_step.back());
    
    Eigen::Vector2cd A_z = a_.inverse() * z;
    
    return (std::pow(k, (3.0 - 2.0*v)) / (2.0 * M_PI * M_PI)) * std::pow(std::abs((std::pow(2, v) * Gamma(v) / (A_z(0) * M_PI)) * (c - d)), 2);

}

Eigen::Matrix2d ModeSolver::A(double x, double p)
{
    Eigen::Matrix2d A;
    double Ai, Bi, Aip, Bip;
    Airy(x, Ai, Aip, Bi, Bip);
    Aip *= p;
    Bip *= p;

    A  << Ai,   Bi,
          Aip,  Bip;

    return A;
}

Eigen::Matrix2cd ModeSolver::H(double eta, double k)
{
    auto x = k*eta;
    auto p = sqrt(eta);
    auto H1 = p * Hankel1(0,x);
    auto H2 = p * Hankel2(0,x);
    auto dH1 = 0.5 * Hankel1(0,x) / p - k * p * Hankel1(1,x);
    auto dH2 = 0.5 * Hankel2(0,x) / p - k * p * Hankel2(1,x);

    Eigen::Matrix2cd H;
    H  << H1,  H2,
          dH1, dH2;

    return H;
}

Eigen::Matrix2cd ModeSolver::H(double x, double k, double v)
{
    auto H1 = I * sqrt(-x) * Hankel1(v, k * x);
    auto H2 = I * sqrt(-x) * Hankel2(v, k * x);
    auto dH1 = 0.5 * I * (sqrt(-x) * k * (Hankel1(v-1, k * x) - Hankel1(v+1, k * x)) - Hankel1(v, k * x) / sqrt(-x));
    auto dH2 = 0.5 * I * (sqrt(-x) * k * (Hankel2(v-1, k * x) - Hankel2(v+1, k * x)) - Hankel2(v, k * x) / sqrt(-x));

    Eigen::Matrix2cd H;
    H  << H1,  H2,
          dH1, dH2;

    return H;
}

Eigen::Matrix2cd ModeSolver::a(std::complex<double> x, double v)
{
    auto a11 = std::pow(x, (0.5 - v));
    auto a12 = std::pow(x, (0.5 + v));
    auto a21 = (0.5 - v) * std::pow(x, (-0.5 - v));
    auto a22 = (0.5 + v) * std::pow(x, (-0.5 + v));

    Eigen::Matrix2cd a;
    a  << a11, a12,
          a21, a22;

    return a;
}
