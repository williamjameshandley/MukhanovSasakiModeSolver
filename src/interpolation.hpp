#pragma once
#include<map>
#include<cmath>
#include<iostream>

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

template<typename X, typename Y>
struct SemiLogInterpolator : public std::map<X,std::pair<Y,int>>
{
    typedef std::map<X,std::pair<Y,int>> map;

    X max() const { return map::rbegin()->first;}
    X min() const { return map::begin()->first;}

    bool insert(typename map::iterator iter, X x, Y y_, Y lim)
    {
        auto x0 = iter->first;
        auto y0 = iter->second.first;
        auto x1 = std::next(iter)->first;
        auto y1 = std::next(iter)->second.first;

        auto y = ( y0 * (x1-x) + y1 * (x-x0) ) / (x1-x0); 
        auto err = std::abs((y-y_)/y_);

        if (y0>0 and y1>0)
        {
            auto yp =  std::exp((std::log(y0) * (x1-x) + std::log(y1) * (x-x0) ) / (x1-x0));
            auto errp = std::abs((yp-y_)/y_);
            if (err < lim or errp < lim)
            {
                iter->second.second =  err < errp ? 0 : 1; 
                return true;
            }
        }
        else if (y0<0 and y1<0)
        {
            auto yn = -std::exp((std::log(-y0) * (x1-x) + std::log(-y1) * (x-x0) ) / (x1-x0));
            auto errp = std::abs((yn-y_)/y_);
            if (err < lim or errp < lim)
            {
                iter->second.second =  err < errp ? 0 : -1; 
                return true;
            }
        }
        else if (err < lim)
        {
            iter->second.second = 0;
            return true;
        }
        map::insert(iter,std::make_pair(x,std::make_pair(y_,NAN)));
        return false;
    }

    Y operator() (X x) const
    { 
        if (map::size()==0) return std::nan("");

        typename map::const_iterator search;
        if (x>=max())         search = --map::end();
        else if (x<=min())    search = ++map::begin();
        else                   search = map::lower_bound(x);

        auto x1 = search->first;
        auto y1 = search->second.first;
        --search;
        auto x0 = search->first;
        auto y0 = search->second.first;
        auto i = search->second.second;

        if (i==0)
            return ( y0 * (x1-x) + y1 * (x-x0) ) / (x1-x0);
        else if (i==1 or i==-1)
            return i*std::exp((std::log(i*y0) * (x1-x) + std::log(i*y1) * (x-x0) ) / (x1-x0));
        else
            throw std::runtime_error("SemiLogInterpolator: type not set");
    }

};
