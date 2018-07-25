#include <iostream>
#include <fstream>
#include "src/BackgroundSolver.cpp"
#include "src/ModeSolver.cpp"
#include "src/Transitions.cpp"
#include "src/Potential.hpp"
#include "src/Special_Functions.hpp"
#include "src/linear_interpolation.hpp"
#include <utility>

using RKCP54 = boost::numeric::odeint::runge_kutta_cash_karp54<std::vector<double>>;

int main()
{
    std::cout<<"Solving for Background"<<std::endl;
    //Background Solver Constructor
    BackgroundSolver solver;
    
    //Background Initial Conditions
    double t0 = 6.48757e-6, t1 = 20.0, phi_p = 23.08546, dphi_p = -sqrt(2.0/3.0) / t0;
    
    //Set Potential
    Polynomial pot(1.0);
    
    //Set Integrator
    double abs_err = 1.0e-5, rel_err = 1e-3;
    boost::numeric::odeint::controlled_runge_kutta<RKCP54> integrator(abs_err, rel_err);
    
    //Solve Background Variables
    auto background_sols = solver.Solve(integrator, pot, t0, t1, phi_p, dphi_p);
    
    //////////////////////////////////////////////////////////////////////////////////
    std::cout<<"Finding Transitions: ";
    //Transitions Conditions
    double eta_end = background_sols.eta.back();
    double eta_i = 0.02 * eta_end;
    double eta_f = 0.95 * eta_end;
    
    //Transitions Constructor
    Transitions T(eta_i, eta_f, background_sols);
    
    //Find Transitions
    double error = 1e-4;
    auto transitions_sols = T.Find(error);
    
    std::cout<<transitions_sols.eta_step.size()<<std::endl;
    
    //////////////////////////////////////////////////////////////////////////////////
    std::cout<<"Finding PPS: ";
    
    ModeSolver ms(background_sols, transitions_sols);
    ms.Initial_Conditions(BD, 0.1 * eta_end);
    
    std::vector<std::pair<double, double>> k_pair;
    double k0 = 1.0, k1 = 1.0e6;
    k_pair.push_back(std::make_pair(k0, k1));
    
    LinearInterpolator<double, double> PS;
    PS.insert(k_pair[0].first, ms.PPS(k_pair[0].first));
    PS.insert(k_pair[0].second, ms.PPS(k_pair[0].second));
    
    auto count = 0;
    
    while(k_pair.size() != 0)
    {
        for(size_t n = 0; n < k_pair.size(); n++)
        {
            auto k_m = exp((log(k_pair[n].first) + log(k_pair[n].second)) / 2.0);
            auto temp_true = ms.PPS(k_m);
            count += 1;
            auto temp_approx = PS(k_m);
            k_pair.erase(k_pair.begin() + n);
            if(abs(temp_true - temp_approx) / temp_true > 0.0005)
            {
                k_pair.push_back(std::make_pair(k0, k_m));
                k_pair.push_back(std::make_pair(k_m, k1));
            }
            PS.insert(k_m, temp_true);
        }
    }
    std::cout<<count<<std::endl;
    
    //////////////////////////////////////////////////////////////////////////////////
    std::cout<<"Plotting..."<<std::endl;
    size_t N = 10000;
    std::vector<double> kplot(N);
    
    for(size_t n = 0; n < N; n++)
    {
        kplot[n] = exp(n * 1.0 * (log(k1) - log(k0)) / N);
    }
    
    std::ofstream mout;
    mout.open("output/PPS.txt");
    for(size_t n = 0; n < N; n++)
    {
        mout<<kplot[n]<<"   "<<PS(kplot[n])<<std::endl;
    }
    mout.close();
    
    
    return 0;
}
