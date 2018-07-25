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
    double k0 = 1.0, k1 = 1.0e6, lim = 1e-3;
    k_pair.push_back(std::make_pair(k0, k1));
    
    LinearInterpolator<double, double> PS;
    PS.insert(k_pair[0].first, ms.PPS(k_pair[0].first));
    PS.insert(k_pair[0].second, ms.PPS(k_pair[0].second));
    
    auto count = 0;
    std::vector<double> kplot;
    while(k_pair.size() != 0)
    {
        for(size_t n = 0; n < k_pair.size(); n++)
        {
            count += 1;
            k0 = k_pair[n].first;
            k1 = k_pair[n].second;
            auto k_m1 = exp((2 * log(k0) + log(k1)) / 3.0);
            auto k_m2 = exp((log(k0) + 2 * log(k1)) / 3.0);
            kplot.push_back(k_m1);
            kplot.push_back(k_m2);
            auto temp_true1 = log(ms.PPS(k_m1));
            auto temp_approx1 = log(PS(k_m1));
            auto temp_true2 = log(ms.PPS(k_m2));
            auto temp_approx2 = log(PS(k_m2));
            k_pair.erase(k_pair.begin() + n);
            if(abs(temp_true1 - temp_approx1) / temp_true1 > lim and abs(temp_true2 - temp_approx2) / temp_true2 > lim)
            {
                k_pair.insert(k_pair.begin() + n, std::make_pair(k0, k_m1));
                k_pair.insert(k_pair.begin() + n + 1, std::make_pair(k_m1, k_m2));
                k_pair.insert(k_pair.begin() + n + 2, std::make_pair(k_m2, k1));
                n += 2;
            }
            else if(abs(temp_true1 - temp_approx1) / temp_true1 > lim and abs(temp_true2 - temp_approx2) / temp_true2 < lim)
            {
                k_pair.insert(k_pair.begin() + n, std::make_pair(k0, k_m1));
                k_pair.insert(k_pair.begin() + n + 1, std::make_pair(k_m1, k1));
                n += 1;
            }
            else if(abs(temp_true1 - temp_approx1) / temp_true1 < lim and abs(temp_true2 - temp_approx2) / temp_true2 > lim)
            {
                k_pair.insert(k_pair.begin() + n, std::make_pair(k0, k_m2));
                k_pair.insert(k_pair.begin() + n + 1, std::make_pair(k_m2, k1));
                n += 1;
            }
            PS.insert(k_m1, exp(temp_true1));
            PS.insert(k_m2, exp(temp_true2));
        }
    }
    std::cout<<count<<std::endl;
    
    //////////////////////////////////////////////////////////////////////////////////
    std::cout<<"Plotting..."<<std::endl;
    
    std::sort(kplot.begin(), kplot.end());
    
    std::ofstream mout;
    mout.open("output/PPS.txt");
    for(size_t n = 0; n < kplot.size(); n++)
    {
        mout<<kplot[n]<<"   "<<PS(kplot[n])<<std::endl;
    }
    mout.close();
    
    
    return 0;
}
