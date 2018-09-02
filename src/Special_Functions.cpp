#include "Special_Functions.hpp"
#include <iostream>

std::complex<double> Bessel_J(int v, double x)
{
    if(x >= 0) return jv(v,x);
    else       return pow(-1, v) * jv(v,-x);
}

std::complex<double> Bessel_Y(int v, double x)
{
    if(x >= 0) return yv(v,x);
    else       return pow(-1, v) * yv(v,-x) + I * 2.0 * pow(-1, v) * Bessel_J(v, -x);
}

std::complex<double> Bessel_I(int v, double x)
{
    if(x >= 0) return iv(v,x);
    else       return pow(-1, v) * iv(v,-x);
}

std::complex<double> Bessel_K(int v, double x)
{
    if(x >= 0) return kn(v,x);
    else       return pow(-1, v) * kn(v,-x) - I * M_PI * Bessel_I(v, -x);
}

void Airy(double x, double& Ai, double& Aip, double& Bi, double& Bip)
{
    airy(x,&Ai,&Aip,&Bi,&Bip);
}

