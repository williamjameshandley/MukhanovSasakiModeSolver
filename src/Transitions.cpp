#include "Transitions.hpp"

double Transitions::True(double N, double k)
{
    return Bsol.omega_2(N) + k * k / exp(2 * Bsol.log_aH(N));
}

double Transitions::Find_N_log_end(double k, double N_i)
{
    std::vector<double> N(10000);
    for(size_t o = 0; o < N.size(); o++)
    {
        N[o] = 1.0 * o * Bsol.N_end / N.size();
    }
    //Find Log to Lin start point
    double N_log_end = N_i;
    if(True(N[0], k) > 0)
    {
        double old = 0;
        int count = 0;
        for(size_t o = 0; o < N.size(); o++)
        {
            double New = True(N[o], k) - 1;
            if(New < 0 and old > 0 and count == 0 and N[o-1] > N_i)
            {
                N_log_end = N[o-1];
                count += 1;
            }
            old = New;
        }
    }
    return N_log_end;
}

TransitionsSolution Transitions::Find(double k, double error)
{
    double N_log_end = Find_N_log_end(k, N_initial);
    
    //Find Log
    std::vector<double> log_N_step;
    std::vector<double> log_a, log_b;
    if(N_log_end != N_initial)
    {
        log_N_step = Log(k, N_initial, N_log_end, error);
        for(size_t l = 0; l < log_N_step.size() - 1; l++)
        {
            log_b.push_back(0.5 * (log(True(log_N_step[l+1], k)) - log(True(log_N_step[l], k))) / (log_N_step[l+1]-log_N_step[l]));
            log_a.push_back(-log_b[l] * log_N_step[l] + 0.5 * log(True(log_N_step[l], k)));
        }
    }
    
    //Find Lin
    std::vector<double> lin_N_step = Linear(k, N_log_end, N_final, error);
    std::vector<double> lin_a, lin_b;
    for(size_t l = 0; l < lin_N_step.size() - 1; l++)
    {
        lin_b.push_back((True(lin_N_step[l+1], k) - True(lin_N_step[l], k)) / (lin_N_step[l+1]-lin_N_step[l]));
        lin_a.push_back(-lin_b[l] * lin_N_step[l] + True(lin_N_step[l], k));
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

std::vector<double> Transitions::Linear(double k, double N_i, double N_f, double lim)
{
    std::vector<double> N_step;
    std::vector<std::pair<double, double>> N_pair;
    LinearInterpolator<double, double> lin_Seg;
    
    //Initialize Extrema
    if(Bsol.N_extrema.size() != 0)
    {
        auto p0 = static_cast<size_t>(std::lower_bound(Bsol.N_extrema.begin(), Bsol.N_extrema.end(), N_i) - Bsol.N_extrema.begin());
        auto p1 = static_cast<size_t>(std::lower_bound(Bsol.N_extrema.begin(), Bsol.N_extrema.end(), N_f) - Bsol.N_extrema.begin());
        if(p1 <= Bsol.N_extrema.size() and p1 > p0)
        {
            N_pair.push_back(std::make_pair(N_i, Bsol.N_extrema[p0]));
            for(size_t n = p0+1; n < p1; n++)
                N_pair.push_back(std::make_pair(Bsol.N_extrema[n-1], Bsol.N_extrema[n]));
            
            N_pair.push_back(std::make_pair(Bsol.N_extrema[p1-1], N_f));
        }
    }
    
    if(N_pair.size() == 0)
    {
        N_pair.push_back(std::make_pair(N_i, N_f));
    }
    
    for(size_t n = 0; n < N_pair.size(); n++)
    {
        lin_Seg.insert(N_pair[n].first, True(N_pair[n].first, k));
        N_step.push_back(N_pair[n].first);
    }
    lin_Seg.insert(N_pair.back().second, True(N_pair.back().second, k));
    N_step.push_back(N_pair.back().second);
    
    while(N_pair.size() != 0)
    {
        for(size_t n = 0; n < N_pair.size(); n++)
        {
            N_i = N_pair[n].first;
            N_f = N_pair[n].second;
            
            auto N_m1 = (2 * N_i + N_f) / 3.0;
            auto N_m2 = (N_i + 2 * N_f) / 3.0;
            
            auto temp_approx1 = lin_Seg(N_m1);
            auto temp_approx2 = lin_Seg(N_m2);
            auto temp_true1 = True(N_m1, k);
            auto temp_true2 = True(N_m2, k);
            
            N_pair.erase(N_pair.begin() + static_cast<int>(n));
            
            if(abs(temp_true1 - temp_approx1) / abs(temp_true1) > lim and abs(temp_true2 - temp_approx2) / abs(temp_true2) > lim)
            {
                N_pair.insert(N_pair.begin() + static_cast<int>(n), std::make_pair(N_i, N_m1));
                N_pair.insert(N_pair.begin() + static_cast<int>(n) + 1, std::make_pair(N_m1, N_m2));
                N_pair.insert(N_pair.begin() + static_cast<int>(n) + 2, std::make_pair(N_m2, N_f));
                n += 2;
                lin_Seg.insert(N_m1, (temp_true1));
                lin_Seg.insert(N_m2, (temp_true2));
                N_step.push_back(N_m1);
                N_step.push_back(N_m2);
            }
            if(abs(temp_true1 - temp_approx1) / abs(temp_true1) > lim and abs(temp_true2 - temp_approx2) / abs(temp_true2) < lim)
            {
                N_pair.insert(N_pair.begin() + static_cast<int>(n), std::make_pair(N_i, N_m1));
                N_pair.insert(N_pair.begin() + static_cast<int>(n) + 1, std::make_pair(N_m1, N_m2));
                n += 1;
                lin_Seg.insert(N_m1, (temp_true1));
                N_step.push_back(N_m1);
            }
            if(abs(temp_true1 - temp_approx1) / abs(temp_true1) < lim and abs(temp_true2 - temp_approx2) / abs(temp_true2) > lim)
            {
                N_pair.insert(N_pair.begin() + static_cast<int>(n), std::make_pair(N_m1, N_m2));
                N_pair.insert(N_pair.begin() + static_cast<int>(n) + 1, std::make_pair(N_m2, N_f));
                n += 1;
                lin_Seg.insert(N_m2, (temp_true2));
                N_step.push_back(N_m2);
            }
        }
    }
    std::sort(N_step.begin(), N_step.end());
    
    return N_step;
    
}

std::vector<double> Transitions::Log(double k, double N_i, double N_f, double lim)
{
    std::vector<double> N_step;
    std::vector<std::pair<double, double>> N_pair;
    LinearInterpolator<double, double> log_Seg;
    
    //Initialize Extrema
    if(Bsol.N_extrema.size() != 0)
    {
        auto p0 = static_cast<size_t>(std::lower_bound(Bsol.N_extrema.begin(), Bsol.N_extrema.end(), N_i) - Bsol.N_extrema.begin());
        auto p1 = static_cast<size_t>(std::lower_bound(Bsol.N_extrema.begin(), Bsol.N_extrema.end(), N_f) - Bsol.N_extrema.begin());
        if(p1 <= Bsol.N_extrema.size() and p1 > p0)
        {
            N_pair.push_back(std::make_pair(N_i, Bsol.N_extrema[p0]));
            for(size_t n = p0+1; n < p1; n++)
                N_pair.push_back(std::make_pair(Bsol.N_extrema[n-1], Bsol.N_extrema[n]));
            
            N_pair.push_back(std::make_pair(Bsol.N_extrema[p1-1], N_f));
        }
    }
    
    if(N_pair.size() == 0)
    {
        N_pair.push_back(std::make_pair(N_i, N_f));
    }
    
    for(size_t n = 0; n < N_pair.size(); n++)
    {
        log_Seg.insert(N_pair[n].first, log(True(N_pair[n].first, k)));
        N_step.push_back(N_pair[n].first);
    }
    log_Seg.insert(N_pair.back().second, log(True(N_pair.back().second, k)));
    N_step.push_back(N_pair.back().second);
    
    while(N_pair.size() != 0)
    {
        for(size_t n = 0; n < N_pair.size(); n++)
        {
            N_i = N_pair[n].first;
            N_f = N_pair[n].second;
            
            auto N_m1 = (2 * N_i + N_f) / 3.0;
            auto N_m2 = (N_i + 2 * N_f) / 3.0;
            
            auto temp_approx1 = log_Seg(N_m1);
            auto temp_approx2 = log_Seg(N_m2);
            auto temp_true1 = log(True(N_m1, k));
            auto temp_true2 = log(True(N_m2, k));
            
            N_pair.erase(N_pair.begin() + static_cast<int>(n));
            
            if(abs(temp_true1 - temp_approx1) / abs(temp_true1) > lim and abs(temp_true2 - temp_approx2) / abs(temp_true2) > lim)
            {
                N_pair.insert(N_pair.begin() + static_cast<int>(n), std::make_pair(N_i, N_m1));
                N_pair.insert(N_pair.begin() + static_cast<int>(n) + 1, std::make_pair(N_m1, N_m2));
                N_pair.insert(N_pair.begin() + static_cast<int>(n) + 2, std::make_pair(N_m2, N_f));
                n += 2;
                log_Seg.insert(N_m1, (temp_true1));
                log_Seg.insert(N_m2, (temp_true2));
                N_step.push_back(N_m1);
                N_step.push_back(N_m2);
            }
            if(abs(temp_true1 - temp_approx1) / abs(temp_true1) > lim and abs(temp_true2 - temp_approx2) / abs(temp_true2) < lim)
            {
                N_pair.insert(N_pair.begin() + static_cast<int>(n), std::make_pair(N_i, N_m1));
                N_pair.insert(N_pair.begin() + static_cast<int>(n) + 1, std::make_pair(N_m1, N_m2));
                n += 1;
                log_Seg.insert(N_m1, (temp_true1));
                N_step.push_back(N_m1);
            }
            if(abs(temp_true1 - temp_approx1) / abs(temp_true1) < lim and abs(temp_true2 - temp_approx2) / abs(temp_true2) > lim)
            {
                N_pair.insert(N_pair.begin() + static_cast<int>(n), std::make_pair(N_m1, N_m2));
                N_pair.insert(N_pair.begin() + static_cast<int>(n) + 1, std::make_pair(N_m2, N_f));
                n += 1;
                log_Seg.insert(N_m2, (temp_true2));
                N_step.push_back(N_m2);
            }
        }
    }
    std::sort(N_step.begin(), N_step.end());
    
    return N_step;
    
}

