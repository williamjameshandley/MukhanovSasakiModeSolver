#include "Transitions.hpp"

TransitionsSolution Transitions::Find(double k, double error)
{
    //Interpolate omega_2
    LinearInterpolator<double, double> True;
    for(size_t o = 0; o < Bsol.N.size(); o++)
        True.insert(Bsol.N[o], Bsol.omega_2[o] + pow(k * exp(-Bsol.N[o]) / Bsol.H[o], 2));
    
    ////////////////////////////////FIX THISSSSSSS
    //Find Log to Lin start point
    double N_log_end = N_i;
    double old = 0;
    int count = 0;
    for(size_t o = 0; o < Bsol.N.size(); o++)
    {
        double New = True(Bsol.N[o]) - 1;
        if(New < 0 and old > 0 and Bsol.d_omega_2[o] < 0 and True(Bsol.N[0]) > 0 and count == 0 and Bsol.N[o] > N_i and True(Bsol.N[o]) > 0)
        {
            N_log_end = Bsol.N[o];
            count += 1;
        }
        old = New;
    }
    if(N_log_end > 10)
    {
        N_log_end = 10;
    }
    
    //Find Log
    std::vector<double> log_N_step;
    std::vector<double> log_a, log_b;
    if(N_log_end != N_i)
    {
        log_N_step = Log(N_i, N_log_end, True, error);
        for(size_t l = 0; l < log_N_step.size() - 1; l++)
        {
            log_b.push_back((log(True(log_N_step[l+1])) - log(True(log_N_step[l]))) / (log_N_step[l+1]-log_N_step[l]));
            log_a.push_back(-log_b[l] * log_N_step[l] + log(True(log_N_step[l])));
        }
    }
    
    //Find Lin
    std::vector<double> lin_N_step = Linear(N_log_end, N_f, True, error);
    std::vector<double> lin_a, lin_b;
    for(size_t l = 0; l < lin_N_step.size() - 1; l++)
    {
        lin_b.push_back((True(lin_N_step[l+1]) - True(lin_N_step[l])) / (lin_N_step[l+1]-lin_N_step[l]));
        lin_a.push_back(-lin_b[l] * lin_N_step[l] + True(lin_N_step[l]));
    }
    
    TransitionsSolution Tsol;
    Tsol.lin_a = lin_a;
    Tsol.lin_b = lin_b;
    Tsol.lin_N_step = lin_N_step;
    Tsol.log_a = log_a;
    Tsol.log_b = log_b;
    Tsol.log_N_step = log_N_step;

    
    return Tsol;
    
}

std::vector<double> Transitions::Linear(double N_initial, double N_final, LinearInterpolator<double, double> True, double lim)
{
    std::vector<double> N_step;
    std::vector<std::pair<double, double>> N_pair;
    LinearInterpolator<double, double> Seg;
    
    auto p0 = static_cast<size_t>(std::lower_bound(Bsol.N.begin(), Bsol.N.end(), N_initial) - Bsol.N.begin());
    auto p1 = static_cast<size_t>(std::lower_bound(Bsol.N.begin(), Bsol.N.end(), N_final) - Bsol.N.begin());
    
    double N0 = N_initial;
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
        N_pair.push_back(std::make_pair(N_initial, N_final));
    }
    else
    {
        N_pair.push_back(std::make_pair(Bsol.N[pp], N_final));
    }
    
    for(size_t n = 0; n < N_pair.size(); n++)
    {
        Seg.insert(N_pair[n].first, True(N_pair[n].first));
        N_step.push_back(N_pair[n].first);
    }
    
    Seg.insert(N_pair.back().second, True(N_pair.back().second));
    N_step.push_back(N_pair.back().second);
    
    while(N_pair.size() != 0)
    {
        for(size_t n = 0; n < N_pair.size(); n++)
        {
            N_initial = N_pair[n].first;
            N_final = N_pair[n].second;
            
            auto N_m1 = (2 * N_initial + N_final) / 3.0;
            auto N_m2 = (N_initial + 2 * N_final) / 3.0;
            
            auto temp_approx1 = (Seg(N_m1));
            auto temp_approx2 = (Seg(N_m2));
            auto temp_true1 = (True(N_m1));
            auto temp_true2 = (True(N_m2));
            
            N_pair.erase(N_pair.begin() + static_cast<int>(n));
            
            if(abs(temp_true1 - temp_approx1) / abs(temp_true1) > lim and abs(temp_true2 - temp_approx2) / abs(temp_true2) > lim)
            {
                N_pair.insert(N_pair.begin() + static_cast<int>(n), std::make_pair(N_initial, N_m1));
                N_pair.insert(N_pair.begin() + static_cast<int>(n) + 1, std::make_pair(N_m1, N_m2));
                N_pair.insert(N_pair.begin() + static_cast<int>(n) + 2, std::make_pair(N_m2, N_final));
                n += 2;
                Seg.insert(N_m1, (temp_true1));
                Seg.insert(N_m2, (temp_true2));
                N_step.push_back(N_m1);
                N_step.push_back(N_m2);
            }
            if(abs(temp_true1 - temp_approx1) / abs(temp_true1) > lim and abs(temp_true2 - temp_approx2) / abs(temp_true2) < lim)
            {
                N_pair.insert(N_pair.begin() + static_cast<int>(n), std::make_pair(N_initial, N_m1));
                N_pair.insert(N_pair.begin() + static_cast<int>(n) + 1, std::make_pair(N_m1, N_final));
                n += 1;
                Seg.insert(N_m1, (temp_true1));
                N_step.push_back(N_m1);
            }
            if(abs(temp_true1 - temp_approx1) / abs(temp_true1) < lim and abs(temp_true2 - temp_approx2) / abs(temp_true2) > lim)
            {
                N_pair.insert(N_pair.begin() + static_cast<int>(n), std::make_pair(N_initial, N_m1));
                N_pair.insert(N_pair.begin() + static_cast<int>(n) + 1, std::make_pair(N_m2, N_final));
                n += 1;
                Seg.insert(N_m2, (temp_true2));
                N_step.push_back(N_m2);
            }
        }
    }
    std::sort(N_step.begin(), N_step.end());
    
    return N_step;
    
}

std::vector<double> Transitions::Log(double N_initial, double N_final, LinearInterpolator<double, double> True, double lim)
{
    std::vector<double> N_step;
    std::vector<std::pair<double, double>> N_pair;
    LinearInterpolator<double, double> Seg;
    
    auto p0 = static_cast<size_t>(std::lower_bound(Bsol.N.begin(), Bsol.N.end(), N_initial) - Bsol.N.begin());
    auto p1 = static_cast<size_t>(std::lower_bound(Bsol.N.begin(), Bsol.N.end(), N_final) - Bsol.N.begin());
    
    double N0 = N_initial;
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
        N_pair.push_back(std::make_pair(N_initial, N_final));
    }
    else
    {
        N_pair.push_back(std::make_pair(Bsol.N[pp], N_final));
    }
    
    for(size_t n = 0; n < N_pair.size(); n++)
    {
        Seg.insert(N_pair[n].first, True(N_pair[n].first));
        N_step.push_back(N_pair[n].first);
    }
    
    Seg.insert(N_pair.back().second, True(N_pair.back().second));
    N_step.push_back(N_pair.back().second);
    
    while(N_pair.size() != 0)
    {
        for(size_t n = 0; n < N_pair.size(); n++)
        {
            N_initial = N_pair[n].first;
            N_final = N_pair[n].second;
            
            auto N_m1 = (2 * N_initial + N_final) / 3.0;
            auto N_m2 = (N_initial + 2 * N_final) / 3.0;
            
            auto temp_approx1 = log(Seg(N_m1));
            auto temp_approx2 = log(Seg(N_m2));
            auto temp_true1 = log(True(N_m1));
            auto temp_true2 = log(True(N_m2));
            
            N_pair.erase(N_pair.begin() + static_cast<int>(n));
            
            if(abs(temp_true1 - temp_approx1) > lim and abs(temp_true2 - temp_approx2) > lim)
            {
                N_pair.insert(N_pair.begin() + static_cast<int>(n), std::make_pair(N_initial, N_m1));
                N_pair.insert(N_pair.begin() + static_cast<int>(n) + 1, std::make_pair(N_m1, N_m2));
                N_pair.insert(N_pair.begin() + static_cast<int>(n) + 2, std::make_pair(N_m2, N_final));
                n += 2;
                Seg.insert(N_m1, exp(temp_true1));
                Seg.insert(N_m2, exp(temp_true2));
                N_step.push_back(N_m1);
                N_step.push_back(N_m2);
            }
            if(abs(temp_true1 - temp_approx1) > lim and abs(temp_true2 - temp_approx2) < lim)
            {
                N_pair.insert(N_pair.begin() + static_cast<int>(n), std::make_pair(N_initial, N_m1));
                N_pair.insert(N_pair.begin() + static_cast<int>(n) + 1, std::make_pair(N_m1, N_final));
                n += 1;
                Seg.insert(N_m1, exp(temp_true1));
                N_step.push_back(N_m1);
            }
            if(abs(temp_true1 - temp_approx1) < lim and abs(temp_true2 - temp_approx2) > lim)
            {
                N_pair.insert(N_pair.begin() + static_cast<int>(n), std::make_pair(N_initial, N_m1));
                N_pair.insert(N_pair.begin() + static_cast<int>(n) + 1, std::make_pair(N_m2, N_final));
                n += 1;
                Seg.insert(N_m2, exp(temp_true2));
                N_step.push_back(N_m2);
            }
        }
    }
    std::sort(N_step.begin(), N_step.end());
    
    return N_step;
    
}

