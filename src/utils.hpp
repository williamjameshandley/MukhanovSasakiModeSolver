#pragma once
#include <functional>
#include <iostream>

template <typename T>
inline T find_root(std::function<T(T)> f, T a, T b, T lim)
{
    auto fa = f(a), fb = f(b);
    if (fa*fb >= 0){
        throw std::runtime_error("find_root: root is not bracketed");
    }

    while ((b-a) > lim)
    {
        auto c = (a+b)/2;
        auto fc = f(c);
        if (fc*fa >= 0) { a=c; fa=fc;}
        else { b=c; fb=fc;} 
    }
    return a;
}
