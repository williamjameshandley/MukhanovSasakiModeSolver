#include "ModeSolver.hpp"
 
void BasicModeSolver::Construct_PPS(double k_i, double k_f, double error = 1e-3)
{
    std::vector<std::pair<double, double>> k_pair;
    
    k_pair.push_back(std::make_pair(k_i, k_f));
    
    PPS.insert(k_pair[0].first, (Find_PPS(k_pair[0].first)));
    PPS.insert(k_pair[0].second, (Find_PPS(k_pair[0].second)));
    
    double lim = error;
    while(k_pair.size() != 0)
    {
        for(size_t n = 0; n < k_pair.size(); n++)
        {
            k_i = k_pair[n].first;
            k_f = k_pair[n].second;
            
            auto k_m1 = exp((3 * log(k_i) + log(k_f)) / 4.0);
            auto k_m2 = exp((log(k_i) + log(k_f)) / 2.0);
            auto k_m3 = exp((log(k_i) + 3 * log(k_f)) / 4.0);
            
            auto temp_approx1 = (PPS(k_m1));
            auto temp_approx2 = (PPS(k_m2));
            auto temp_approx3 = (PPS(k_m3));
            auto temp_true1 = (Find_PPS(k_m1));
            auto temp_true2 = (Find_PPS(k_m2));
            auto temp_true3 = (Find_PPS(k_m3));
            
            k_pair.erase(k_pair.begin() + static_cast<int>(n));
            
            if(abs(temp_true1 - temp_approx1) / abs(temp_true1) > lim and abs(temp_true2 - temp_approx2) / abs(temp_true2) > lim and abs(temp_true3 - temp_approx3) / abs(temp_true3) > lim)
            {
                k_pair.insert(k_pair.begin() + static_cast<int>(n), std::make_pair(k_i, k_m1));
                k_pair.insert(k_pair.begin() + static_cast<int>(n) + 1, std::make_pair(k_m1, k_f));
                k_pair.insert(k_pair.begin() + static_cast<int>(n) + 2, std::make_pair(k_i, k_m2));
                k_pair.insert(k_pair.begin() + static_cast<int>(n) + 3, std::make_pair(k_m2, k_f));
                k_pair.insert(k_pair.begin() + static_cast<int>(n) + 4, std::make_pair(k_i, k_m3));
                k_pair.insert(k_pair.begin() + static_cast<int>(n) + 5, std::make_pair(k_m3, k_f));
                n += 5;
            }
            if(abs(temp_true1 - temp_approx1) / abs(temp_true1) < lim and abs(temp_true2 - temp_approx2) / abs(temp_true2) > lim and abs(temp_true3 - temp_approx3) / abs(temp_true3) > lim)
            {
                k_pair.insert(k_pair.begin() + static_cast<int>(n), std::make_pair(k_i, k_m2));
                k_pair.insert(k_pair.begin() + static_cast<int>(n) + 1, std::make_pair(k_m2, k_f));
                k_pair.insert(k_pair.begin() + static_cast<int>(n) + 2, std::make_pair(k_i, k_m3));
                k_pair.insert(k_pair.begin() + static_cast<int>(n) + 3, std::make_pair(k_m3, k_f));
                n += 3;
            }
            if(abs(temp_true1 - temp_approx1) / abs(temp_true1) > lim and abs(temp_true2 - temp_approx2) / abs(temp_true2) < lim and abs(temp_true3 - temp_approx3) / abs(temp_true3) > lim)
            {
                k_pair.insert(k_pair.begin() + static_cast<int>(n), std::make_pair(k_i, k_m1));
                k_pair.insert(k_pair.begin() + static_cast<int>(n) + 1, std::make_pair(k_m1, k_f));
                k_pair.insert(k_pair.begin() + static_cast<int>(n) + 2, std::make_pair(k_i, k_m3));
                k_pair.insert(k_pair.begin() + static_cast<int>(n) + 3, std::make_pair(k_m3, k_f));
                n += 3;
            }
            if(abs(temp_true1 - temp_approx1) / abs(temp_true1) < lim and abs(temp_true2 - temp_approx2) / abs(temp_true2) < lim and abs(temp_true3 - temp_approx3) / abs(temp_true3) > lim)
            {
                k_pair.insert(k_pair.begin() + static_cast<int>(n), std::make_pair(k_i, k_m3));
                k_pair.insert(k_pair.begin() + static_cast<int>(n) + 1, std::make_pair(k_m3, k_f));
                n += 1;
            }
            if(abs(temp_true1 - temp_approx1) / abs(temp_true1) > lim and abs(temp_true2 - temp_approx2) / abs(temp_true2) > lim and abs(temp_true3 - temp_approx3) / abs(temp_true3) < lim)
            {
                k_pair.insert(k_pair.begin() + static_cast<int>(n), std::make_pair(k_i, k_m1));
                k_pair.insert(k_pair.begin() + static_cast<int>(n) + 1, std::make_pair(k_m1, k_f));
                k_pair.insert(k_pair.begin() + static_cast<int>(n) + 2, std::make_pair(k_i, k_m2));
                k_pair.insert(k_pair.begin() + static_cast<int>(n) + 3, std::make_pair(k_m2, k_f));
                n += 3;
            }
            if(abs(temp_true1 - temp_approx1) / abs(temp_true1) < lim and abs(temp_true2 - temp_approx2) / abs(temp_true2) > lim and abs(temp_true3 - temp_approx3) / abs(temp_true3) < lim)
            {
                k_pair.insert(k_pair.begin() + static_cast<int>(n), std::make_pair(k_i, k_m2));
                k_pair.insert(k_pair.begin() + static_cast<int>(n) + 1, std::make_pair(k_m2, k_f));
                n += 1;
            }
            if(abs(temp_true1 - temp_approx1) / abs(temp_true1) > lim and abs(temp_true2 - temp_approx2) / abs(temp_true2) < lim and abs(temp_true3 - temp_approx3) / abs(temp_true3) < lim)
            {
                k_pair.insert(k_pair.begin() + static_cast<int>(n), std::make_pair(k_i, k_m1));
                k_pair.insert(k_pair.begin() + static_cast<int>(n) + 1, std::make_pair(k_m1, k_f));
                n += 1;
            }
            
            PPS.insert(k_m1, (temp_true1));
            PPS.insert(k_m2, (temp_true2));
            PPS.insert(k_m3, (temp_true3));
            k_plot.push_back(k_m1);
            k_plot.push_back(k_m2);
            k_plot.push_back(k_m3);
        }
    }
    std::sort(k_plot.begin(), k_plot.end());
}
