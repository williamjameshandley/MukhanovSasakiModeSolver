#include "Transitions.hpp"
#include "linear_interpolation.hpp"
#include <boost/math/tools/minima.hpp>
#include <cmath>

//Integral of ddz Function
double Transitions::integral(double a, double b) {
    
    double area = 0.0;
    size_t n_a = static_cast<size_t>(std::lower_bound(eta_sol.begin(), eta_sol.end(), a) - eta_sol.begin());
    size_t n_b = static_cast<size_t>(std::lower_bound(eta_sol.begin(), eta_sol.end(), b) - eta_sol.begin());
    double dx = (b - a) / (n_b - n_a);
    
    area += ddz[n_a] * dx / 2.;
    
    for (size_t i = n_a + 1; i < n_b; i ++) {
        area += ddz[i] * dx;
    }
    
    area += ddz[n_b] * dx / 2.;
    
    return area;
}

std::tuple<std::vector<double>, std::vector<double>, double,  std::vector<double>> Transitions::Find(double error)
{
    double dA = 0, F = 0, aa, bb, delta;
    size_t n = 0;
    std::vector<double> a, b, eta_step;
    
    eta_step.push_back(eta_i);
    size_t i = static_cast<size_t>(std::lower_bound(eta_sol.begin(), eta_sol.end(), eta_i) - eta_sol.begin()) + 1;
    eta_step.push_back(eta_sol[i]);
    
    LinearInterpolator<double, double> DDZ;
    for(size_t o = 0; o < ddz.size(); o++)
    {
        DDZ.insert(eta_sol[o], ddz[o]);
    }
    
    while(eta_step[n+1] < eta_f)
    {
        while(dA < error)
        {
            eta_step[n+1] = eta_sol[i];
            F = integral(eta_step[n], eta_step[n+1]);
            bb = (DDZ(eta_step[n+1]) - DDZ(eta_step[n])) / (eta_step[n+1]-eta_step[n]);
            aa = -bb * eta_step[n] + DDZ(eta_step[n]);
            dA = abs((F - (aa*(eta_step[n+1] - eta_step[n]) + 0.5 * bb * (eta_step[n+1] * eta_step[n+1] - eta_step[n] * eta_step[n]))) / F);
            i += 1;
            
        }
        dA = 0;
        eta_step.push_back(eta_step[n+1]);
        n += 1;
    }
    
    for(size_t l = 0; l < eta_step.size() - 1; l++)
    {
        b.push_back(( DDZ(eta_step[n+1]) -  DDZ(eta_step[n])) / (eta_step[l+1]-eta_step[l]));
        a.push_back(-b[l] * eta_step[l] +  DDZ(eta_step[n]));
    }
    
    delta = (DDZ(eta_step.back()) / (2/((eta_step.back() - eta_end) * (eta_step.back() - eta_end))))  -  1;
    
    return std::make_tuple(a, b, delta, eta_step);
    
}


