#include "Special_Functions.hpp"


std::complex<double> Bessel_J(int v, double x)
{
    if(x > 0)
        return jv(v,x);
    
    else if (x < 0)
        return pow(-1, v) * jv(v,-x);
    
    else
        return 0;
}

std::complex<double> Bessel_Y(int v, double x)
{
    if(x > 0)
        return yv(v,x);
    else if (x < 0)
        return pow(-1, v) * yv(v,-x) + I * 2.0 * pow(-1, v) * Bessel_J(v, -x);
    else
        return 0;
}

std::complex<double> Hankel1(int v, double x)
{
    return Bessel_J(v, x) + I * Bessel_Y(v, x);
}

std::complex<double> Hankel2(int v, double x)
{
    return Bessel_J(v, x) - I * Bessel_Y(v, x);
}

void Airy(double x, double& Ai, double& Aip, double& Bi, double& Bip){airy(x,&Ai,&Aip,&Bi,&Bip);}

