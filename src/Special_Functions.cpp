#include "Special_Functions.hpp"
#include <iostream>


std::complex<double> Bessel_J(double v, double x)
{
    if(x > 0)
        return jv(v,x);
    else if (x < 0)
        return exp(I * v * M_PI) * jv(v,-x);
    else
        return 0;
}

std::complex<double> Bessel_Y(double v, double x)
{
    if(x > 0)
        return yv(v,x);
    else if (x < 0)
        return exp(I * v * M_PI) * yv(v,-x);
    else
        return 0;
}

std::complex<double> Hankel1(double v, double x) { return Bessel_J(v, x) + I * Bessel_Y(v, x); }
std::complex<double> Hankel2(double v, double x) { return Bessel_J(v, x) - I * Bessel_Y(v, x); }

void Airy(double x, double& Ai, double& Aip, double& Bi, double& Bip){airy(x,&Ai,&Aip,&Bi,&Bip);}

double Gamma(double x) { return gamma(x); }
