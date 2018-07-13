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
//    interp.insert(1,1);
//    interp.insert(2,-1);
//    interp.insert(3,1);
//    interp.insert(4,-1);
//    for (auto t=-1.;t<5;t+=0.1)
//        std::cout << t << " " << interp(t) << std::endl;

template<typename X, typename Y>
struct LinearInterpolator
{
    std::map<X,Y> points;

    LinearInterpolator() : points{} {}

    void insert(X x, Y y) {points[x] = y;}

    X xmax() const { return points.rbegin()->first;}
    X xmin() const { return points.begin()->first;}

    Y operator() (X x) const
    { 
        if (points.size()==0) return std::nan("");

        typename std::map<X,Y>::const_iterator search;
        if (x>=xmax())         search = --points.end();
        else if (x<=xmin())    search = ++points.begin();
        else                   search = points.lower_bound(x);

        auto x1 = search->first;
        auto y1 = search->second;
        --search;
        auto x0 = search->first;
        auto y0 = search->second;

        return ( y0 * (x1-x) + y1 * (x-x0) ) / (x1-x0);
    }

};
