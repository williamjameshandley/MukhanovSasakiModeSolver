#include "Transitions.hpp"
#include <algorithm>

double Transitions::True(double N, double k)
{
    return omega_2(N) + k * k / exp(2 * log_aH(N));
}

TransitionsSolution Transitions::Find(double k, double error)
{
    //Find Distribution and Segment control
    auto seg_control = N_Distribution(k, N_initial, N_final, error);
    
    //Find Segments
    std::vector<double> a, b;
    for(size_t n = 0; n < seg_control.size() - 1; n++)
    {
        if(seg_control[n].second == 0)
        {
            b.push_back((True(seg_control[n+1].first, k) - True(seg_control[n].first, k)) / (seg_control[n+1].first - seg_control[n].first));
            a.push_back(True(seg_control[n].first, k) - b[n] * seg_control[n].first);
        }
        else if(seg_control[n].second == 1)
        {
            b.push_back(0.5 * (log(True(seg_control[n+1].first, k)) - log(True(seg_control[n].first, k))) / (seg_control[n+1].first - seg_control[n].first));
            a.push_back(0.5 * log(True(seg_control[n].first, k)) - b[n] * seg_control[n].first);
        }
        else if(seg_control[n].second == -1)
        {
            b.push_back(0.5 * (log(-True(seg_control[n+1].first, k)) - log(-True(seg_control[n].first, k))) / (seg_control[n+1].first - seg_control[n].first));
            a.push_back(0.5 * log(-True(seg_control[n].first, k)) - b[n] * seg_control[n].first);
        }
    }
    
    TransitionsSolution Tsol;
    Tsol.a = a;
    Tsol.b = b;
    Tsol.seg_control = seg_control;
    
    return Tsol;
}

std::vector<std::pair<double, int>> Transitions::N_Distribution(double k, double N_i, double N_f, double lim)
{
    std::vector<std::pair<double, double>> N_pair;
    std::vector<std::pair<double, int>> seg_control;
    seg_control.push_back(std::make_pair(N_f, 0));
    
    //Initialize Extrema
    if(N_extrema.size() != 0)
    {
        auto p0 = static_cast<size_t>(std::lower_bound(N_extrema.begin(), N_extrema.end(), N_i) - N_extrema.begin());
        auto p1 = static_cast<size_t>(std::lower_bound(N_extrema.begin(), N_extrema.end(), N_f) - N_extrema.begin());
        if(p1 <= N_extrema.size() and p1 > p0)
        {
            N_pair.push_back(std::make_pair(N_i, N_extrema[p0]));
            for(size_t n = p0+1; n < p1; n++)
                N_pair.push_back(std::make_pair(N_extrema[n-1], N_extrema[n]));
            
            N_pair.push_back(std::make_pair(N_extrema[p1-1], N_f));
        }
    }
    if(N_pair.size() == 0)
    {
        N_pair.push_back(std::make_pair(N_i, N_f));
    }
    
    while(N_pair.size() != 0)
    {
        for(size_t n = 0; n < N_pair.size(); n++)
        {
            N_i = N_pair[n].first;
            N_f = N_pair[n].second;
            N_pair.erase(N_pair.begin() + static_cast<int>(n));
            auto End0 = True(N_i, k);
            auto End1 = True(N_f, k);
            
            auto N_m = (N_i + N_f) / 2.0;
            auto temp_true = True(N_m, k);
            
            //Try exp or linear
            if(End0 > 1 and End1 > 1)
            {
                auto lin_b = (End1 - End0) / (N_f - N_i);
                auto lin_a = End0 - lin_b * N_i;
                auto lin_approx = lin_a + lin_b * N_m;
                
                auto log_b = (log(End1) - log(End0)) / (N_f - N_i);
                auto log_a = log(End0) - log_b * N_i;
                auto log_approx = log_a + log_b * N_m;
                
                auto err_lin = abs(lin_approx - temp_true) / abs(temp_true);
                auto err_log = abs(exp(log_approx) - temp_true) / abs(temp_true);
                
                if(err_lin > err_log)
                {
                    if(err_log > lim)
                    {
                        N_pair.insert(N_pair.begin() + static_cast<int>(n), std::make_pair(N_i, N_m));
                        N_pair.insert(N_pair.begin() + static_cast<int>(n) + 1, std::make_pair(N_m, N_f));
                        n += 1;
                    }
                    else
                    {
                        seg_control.push_back(std::make_pair(N_i, 1));
                    }
                }
                else if(err_lin < err_log)
                {
                    if(err_lin > lim)
                    {
                        N_pair.insert(N_pair.begin() + static_cast<int>(n), std::make_pair(N_i, N_m));
                        N_pair.insert(N_pair.begin() + static_cast<int>(n) + 1, std::make_pair(N_m, N_f));
                        n += 1;
                    }
                    else
                    {
                        seg_control.push_back(std::make_pair(N_i, 0));
                    }
                }
            }
            //Try -exp or linear
            else if(End0 < -1 and End1 < -1)
            {
                auto lin_b = (End1 - End0) / (N_f - N_i);
                auto lin_a = End0 - lin_b * N_i;
                auto lin_approx = lin_a + lin_b * N_m;

                auto log_b = (log(-End1) - log(-End0)) / (N_f - N_i);
                auto log_a = log(-End0) - log_b * N_i;
                auto log_approx = log_a + log_b * N_m;

                auto err_lin = abs(lin_approx - temp_true) / abs(temp_true);
                auto err_log = abs(-exp(log_approx) - temp_true) / abs(temp_true);

                if(err_lin > err_log)
                {
                    if(err_log > lim)
                    {
                        N_pair.insert(N_pair.begin() + static_cast<int>(n), std::make_pair(N_i, N_m));
                        N_pair.insert(N_pair.begin() + static_cast<int>(n) + 1, std::make_pair(N_m, N_f));
                        n += 1;
                    }
                    else
                    {
                        seg_control.push_back(std::make_pair(N_i, -1));
                    }
                }
                else if(err_lin < err_log)
                {
                    if(err_lin > lim)
                    {
                        N_pair.insert(N_pair.begin() + static_cast<int>(n), std::make_pair(N_i, N_m));
                        N_pair.insert(N_pair.begin() + static_cast<int>(n) + 1, std::make_pair(N_m, N_f));
                        n += 1;
                    }
                    else
                    {
                        seg_control.push_back(std::make_pair(N_i, 0));
                    }
                }
            }
            //linear
            else
            {
                auto lin_b = (End1 - End0) / (N_f - N_i);
                auto lin_a = End0 - lin_b * N_i;
                auto lin_approx = lin_a + lin_b * N_m;

                auto err_lin = abs(lin_approx - temp_true) / abs(temp_true);
                
                if(err_lin > lim)
                {
                    N_pair.insert(N_pair.begin() + static_cast<int>(n), std::make_pair(N_i, N_m));
                    N_pair.insert(N_pair.begin() + static_cast<int>(n) + 1, std::make_pair(N_m, N_f));
                    n += 1;
                }
                else
                {
                    seg_control.push_back(std::make_pair(N_i, 0));
                }
            }
        }
    }
    
    std::sort(seg_control.begin(), seg_control.end());
    return seg_control;
}
