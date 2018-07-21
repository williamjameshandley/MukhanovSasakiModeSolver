#pragma once
#include <iostream>
#include <math.h>
#include <boost/math/special_functions/airy.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>

const std::complex<double> I(0, 1);

inline std::complex<double> Bessel_J(double v, double x)
{
    if(x > 0)
        return boost::math::cyl_bessel_j(v,x);
    else if (x < 0)
        return exp(I * v * M_PI) * boost::math::cyl_bessel_j(v,-x);
    else
        return 0;
}

inline std::complex<double> Bessel_Y(double v, double x)
{
    
    if(x > 0)
        return boost::math::cyl_neumann(v,x);
    else if (x < 0)
        return exp(- I * v * M_PI) * boost::math::cyl_neumann(v,-x) +  I * 2. * cos(M_PI * v) * boost::math::cyl_bessel_j(v,-x);
    else
        return 0;
}

inline std::complex<double> Hankel1(double v, double x) { return Bessel_J(v, x) + I * Bessel_Y(v, x); }
inline std::complex<double> Hankel2(double v, double x) { return Bessel_J(v, x) - I * Bessel_Y(v, x); }

inline double Airy_Ai(double x) { return boost::math::airy_ai(x); }
inline double Airy_Bi(double x) { return boost::math::airy_bi(x); }
inline double Airy_Aip(double x) { return boost::math::airy_ai_prime(x); }
inline double Airy_Bip(double x) { return boost::math::airy_bi_prime(x); }
inline double Gamma(double x) { return tgamma(x); }