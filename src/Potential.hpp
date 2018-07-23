#pragma once
#include <iostream>
#include <math.h>

struct Potential {

    virtual double V(double) =0;
    virtual double dV(double) =0;
    virtual double ddV(double) =0;

};

struct Polynomial : public Potential
{
    Polynomial( double _m, double _lambda=0) : m{_m}, lambda{_lambda} {}

    double m;
    double lambda;

    double V(double phi)   { return 0.5 * m * m * phi * phi + (1.0/24) * lambda * phi * phi * phi * phi; }
    double dV(double phi)  { return m * m * phi + (1.0/6) * lambda * phi * phi * phi; }
    double ddV(double phi) { return m * m + (1.0/2) * lambda * phi * phi; }
};

struct Starobinsky : public Potential
{
    Starobinsky( double _m) : m{_m}  {}

    double m;

    double V(double phi)   { return 0.75 * m * m * (1 - exp(-sqrt(2.0/3) * phi)) * (1 - exp(-sqrt(2.0/3) * phi)); }
    double dV(double phi)  { return 0.75 * m * m * sqrt(8.0/3) * exp(-sqrt(2.0/3) * phi) * (1 - exp(-sqrt(2.0/3) * phi)); }
    double ddV(double phi) { return 0.75 * (4.0/3) * m * m * exp(-2*sqrt(2.0/3) * phi) * (2 - exp(sqrt(2.0/3) * phi)); }
};

