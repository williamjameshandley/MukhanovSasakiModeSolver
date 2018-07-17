#pragma once
#include <iostream>
#include <math.h>
#include <boost/math/special_functions/airy.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/gamma.hpp>

const std::complex<double> ii(0, 1);

std::complex<double> Bessel_J(double v, double x)
{
    std::complex<double> BJ;
    
    if(x > 0)
    {
        BJ = boost::math::cyl_bessel_j(v,x);
    }
    else if (x < 0)
    {
        BJ = exp(ii * v * M_PI) * boost::math::cyl_bessel_j(v,-x);
    }
    
    return BJ;
}

std::complex<double> Bessel_Y(double v, double x)
{
    std::complex<double> BY;
    
    if(x > 0)
    {
        BY = boost::math::cyl_neumann(v,x);
    }
    else if (x < 0)
    {
        BY = exp(- ii * v * M_PI) * boost::math::cyl_neumann(v,-x) +  ii * 2. * cos(M_PI * v) * boost::math::cyl_bessel_j(v,-x);
    }
    
    return BY;
}

double Airy_Ai(double x)
{
    return boost::math::airy_ai(x);
}

double Airy_Bi(double x)
{
    return boost::math::airy_bi(x);
}

double Airy_Aip(double x)
{
    return boost::math::airy_ai_prime(x);
}

double Airy_Bip(double x)
{
    return boost::math::airy_bi_prime(x);
}

double Gamma(double x)
{
    return tgamma(x);
}
