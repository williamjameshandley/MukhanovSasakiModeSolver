#include "Transitions.hpp"
#include "linear_interpolation.hpp"
#include <boost/math/tools/minima.hpp>
#include <cmath>

//Integral of ddz Function
double Transitions::integral(double a, double b) {
    
    auto n_a = static_cast<size_t>(std::lower_bound(Bsol.eta.begin(), Bsol.eta.end(), a) - Bsol.eta.begin());
    auto n_b = static_cast<size_t>(std::lower_bound(Bsol.eta.begin(), Bsol.eta.end(), b) - Bsol.eta.begin());
    auto dx = (b - a) / static_cast<double>(n_b - n_a);
    
    auto area = Bsol.ddz[n_a] * dx / 2.;
    
    for (size_t i = n_a + 1; i < n_b; i ++) 
        area += Bsol.ddz[i] * dx;
    
    area += Bsol.ddz[n_b] * dx / 2.;
    
    return area;
}

TransitionsSolution Transitions::Find(double error)
{
    std::vector<double> eta_step;
    
    eta_step.push_back(eta_i);
    
    size_t p = static_cast<size_t>(std::lower_bound(Bsol.eta.begin(), Bsol.eta.end(), eta_i) - Bsol.eta.begin()) + 1;
    
    eta_step.push_back(Bsol.eta[p]);
    
    LinearInterpolator<double, double> DDZ;
    for(size_t o = 0; o < Bsol.ddz.size(); o++)
        DDZ.insert(Bsol.eta[o], Bsol.ddz[o]);
    
    size_t n = 0;
    while(eta_step[n+1] < eta_f)
    {

        double dA = 0;
        while(dA < error)
        {
            eta_step[n+1] = Bsol.eta[p];
            auto F = integral(eta_step[n], eta_step[n+1]);
            auto bb = (DDZ(eta_step[n+1]) - DDZ(eta_step[n])) / (eta_step[n+1]-eta_step[n]);
            auto aa = -bb * eta_step[n] + DDZ(eta_step[n]);
            dA = abs((F - (aa*(eta_step[n+1] - eta_step[n]) + 0.5 * bb * (eta_step[n+1] * eta_step[n+1] - eta_step[n] * eta_step[n]))) / F);
            p += 1;
            
        }
        
        eta_step.push_back(Bsol.eta[p]);
        n += 1;
    }
    
    std::vector<double> a, b;
    for(size_t l = 0; l < eta_step.size() - 1; l++)
    {
        b.push_back((DDZ(eta_step[l+1]) -  DDZ(eta_step[l])) / (eta_step[l+1]-eta_step[l]));
        a.push_back(-b[l] * eta_step[l] +  DDZ(eta_step[l]));
    }
    
    auto delta = (DDZ(eta_step.back()) / (2/((eta_step.back() - Bsol.eta.back()) * (eta_step.back() - Bsol.eta.back()))))  -  1;
    
    return TransitionsSolution(delta, a, b, eta_step);
    
}


