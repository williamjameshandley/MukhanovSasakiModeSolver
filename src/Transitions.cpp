#include <cmath>
#include <boost/math/tools/minima.hpp>
#include "Transitions.hpp"
#include "linear_interpolation.hpp"

TransitionsSolution Transitions::Find(double k, double error)
{
    LinearInterpolator<double, double> Seg, True;
    std::vector<std::pair<double, double>> N_pair;
    std::vector<double> N_step;
    
    //Interpolate omega_2
    for(size_t o = 0; o < Bsol.N.size(); o++)
        True.insert(Bsol.N[o], Bsol.omega_2[o] + pow(k * exp(-Bsol.N[o]) / Bsol.H[o], 2));
    
    //Initialize pairs at zeros of derivative of omega_2
    auto p0 = static_cast<size_t>(std::lower_bound(Bsol.N.begin(), Bsol.N.end(), N_i) - Bsol.N.begin());
    auto p1 = static_cast<size_t>(std::lower_bound(Bsol.N.begin(), Bsol.N.end(), N_f) - Bsol.N.begin());
    
    double N0 = N_i;
    size_t pp = p0;
    for(size_t o = p0 + 1; o < p1; o++)
    {
        double old = Bsol.d_omega_2[o - 1];
        if(old < 0 and Bsol.d_omega_2[o] > 0)
        {
            N_pair.push_back(std::make_pair(N0, Bsol.N[o]));
            N0 =  Bsol.N[o];
            pp = o;
        }
        else if(old > 0 and Bsol.d_omega_2[o] < 0)
        {
            N_pair.push_back(std::make_pair(N0, Bsol.N[o]));
            N0 =  Bsol.N[o];
            pp = o;
        }
    }
    if(pp == p0)
    {
        N_pair.push_back(std::make_pair(N_i, N_f));
    }
    else
    {
        N_pair.push_back(std::make_pair(Bsol.N[pp], N_f));
    }

    for(size_t n = 0; n < N_pair.size(); n++)
    {
        Seg.insert(N_pair[n].first, True(N_pair[n].first));
        N_step.push_back(N_pair[n].first);
    }
    
    Seg.insert(N_pair.back().second, True(N_pair.back().second));
    N_step.push_back(N_pair.back().second);
    
    double lim = error;
        
    while(N_pair.size() != 0)
    {
        for(size_t n = 0; n < N_pair.size(); n++)
        {
            N_i = N_pair[n].first;
            N_f = N_pair[n].second;
            
            auto N_m1 = (2 * N_i + N_f) / 3.0;
            auto N_m2 = (N_i + 2 * N_f) / 3.0;
            
            auto temp_approx1 = (Seg(N_m1));
            auto temp_approx2 = (Seg(N_m2));
            auto temp_true1 = (True(N_m1));
            auto temp_true2 = (True(N_m2));
            
            N_pair.erase(N_pair.begin() + static_cast<int>(n));
            
            if(abs(temp_true1 - temp_approx1) / abs(temp_true1) > lim or abs(temp_true2 - temp_approx2) / abs(temp_true2) > lim)
            {
                N_pair.insert(N_pair.begin() + static_cast<int>(n), std::make_pair(N_i, N_m1));
                N_pair.insert(N_pair.begin() + static_cast<int>(n) + 1, std::make_pair(N_m1, N_m2));
                N_pair.insert(N_pair.begin() + static_cast<int>(n) + 2, std::make_pair(N_m2, N_f));
                n += 2;
                Seg.insert(N_m1, (temp_true1));
                Seg.insert(N_m2, (temp_true2));
                N_step.push_back(N_m1);
                N_step.push_back(N_m2);
            }
        }
    }
    
    std::sort (N_step.begin(), N_step.end());
    std::vector<double> a, b;
    for(size_t l = 0; l < N_step.size() - 1; l++)
    {
        b.push_back((True(N_step[l+1]) - True(N_step[l])) / (N_step[l+1]-N_step[l]));
        a.push_back(-b[l] * N_step[l] + True(N_step[l]));
    }
    
    
    TransitionsSolution Tsol;
    Tsol.a = a;
    Tsol.b = b;
    Tsol.N_step = N_step;
    
    return Tsol;
    
}


