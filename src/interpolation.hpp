#pragma once
#include<map>
#include<cmath>

// Linear interpolation from type X to type Y
// Works with any vector space Y over scalars X 
// https://en.wikipedia.org/wiki/Vector_space#Definition
//
// example usage:
//
//    LinearInterpolator<double, double> interp;
//    interp[1] = 1;
//    interp[2] = -1;
//    interp[3] = 1;
//    interp[4] = -1;
//    for (auto t=-1.;t<5;t+=0.1)
//        std::cout << t << " " << interp(t) << std::endl;
//
// Derives from std::map.

template<typename X, typename Y>
struct LinearInterpolator : public std::map<X,Y>
{

    X max() const { return std::map<X,Y>::rbegin()->first;}
    X min() const { return std::map<X,Y>::begin()->first;}

    Y operator() (X x) const
    { 
        if (std::map<X,Y>::size()==0) return std::nan("");

        typename std::map<X,Y>::const_iterator search;
        if (x>=max())         search = --std::map<X,Y>::end();
        else if (x<=min())    search = ++std::map<X,Y>::begin();
        else                   search = std::map<X,Y>::lower_bound(x);

        auto x1 = search->first;
        auto y1 = search->second;
        --search;
        auto x0 = search->first;
        auto y0 = search->second;

        return ( y0 * (x1-x) + y1 * (x-x0) ) / (x1-x0);
    }

};
