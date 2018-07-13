#include "Transitions.hpp"
#include <boost/math/tools/minima.hpp>
#include <cmath>

//Integral of ddz Function
double Transitions::integral(double a, double b) {
    
    double area = 0.0;
    size_t n_a = search(eta_sol, a);
    size_t n_b = search(eta_sol, b);
    double dx = (b - a) / (n_b - n_a);
    
    area += ddz[n_a] * dx / 2.;
    
    for (size_t i = n_a + 1; i < n_b; i ++) {
        area += ddz[i] * dx;
    }
    
    area += ddz[n_b] * dx / 2.;
    
    return area;
}

//Vector Search
size_t Transitions::search(std::vector<double> a, double b)
{
    size_t n = 0;
    while (a[n] < b) {n += 1;}
    return n;
}

 std::tuple<std::vector<double>, std::vector<double>, double,  std::vector<double>> Transitions::Find()
{
    double dA = 0, F = 0, aa, bb, delta;
    size_t n = 0;
    std::vector<double> a, b, eta_step;
    
    eta_step.push_back(eta_i);
    size_t i = search(eta_sol, eta_i) + 1;
    eta_step.push_back(eta_sol[i]);
    
    while(eta_step[n+1] < eta_f)
    {
        while(dA < 0.0005)
        {
            eta_step[n+1] = eta_sol[i];
            F = integral(eta_step[n], eta_step[n+1]);
            bb = (ddz[search(eta_sol, eta_step[n+1])] - ddz[search(eta_sol, eta_step[n])]) / (eta_step[n+1]-eta_step[n]);
            aa = -bb * eta_step[n] + ddz[search(eta_sol, eta_step[n])];
            dA = abs((F - (aa*(eta_step[n+1] - eta_step[n]) + 0.5 * bb * (eta_step[n+1] * eta_step[n+1] - eta_step[n] * eta_step[n]))) / F);
            i += 1;
            
        }
        dA = 0;
        eta_step.push_back(eta_step[n+1]);
        n += 1;
        std::cout<<n<<std::endl;
    }
    
    for(size_t l = 0; l < eta_step.size() - 1; l++)
    {
        b.push_back((ddz[search(eta_sol, eta_step[l+1])] - ddz[search(eta_sol, eta_step[l])]) / (eta_step[l+1]-eta_step[l]));
        a.push_back(-b[l] * eta_step[l] + ddz[search(eta_sol, eta_step[l])]);
    }
    
    delta = (ddz[search(eta_sol, eta_step.back())] / (2/((eta_step.back() - eta_end) * (eta_step.back() - eta_end))))  -  1;
    
    return std::make_tuple(a, b, delta, eta_step);
    
}


