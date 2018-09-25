#pragma once
#include <iostream>
#include <math.h>
#include <cmath>

struct Potential {
    double dphi = 1e-8;
    
    virtual double V(double) =0;
    virtual double dV(double phi) =0;
    virtual double ddV(double phi) =0;
    virtual double dddV(double phi) {return (ddV(phi+dphi) - ddV(phi-dphi))/2/dphi;};
    virtual ~Potential() {}

};

struct Polynomial : public Potential
{
    Polynomial( double _m, double _lambda=0) : m{_m}, lambda{_lambda} {}

    double m;
    double lambda;

    double V(double phi)   { return 0.5 * m * m * phi * phi + (1.0/24) * lambda * phi * phi * phi * phi; }
    double dV(double phi)  { return m * m * phi + (1.0/6) * lambda * phi * phi * phi; }
    double ddV(double phi) { return m * m + (1.0/2) * lambda * phi * phi; }
    double dddV(double phi) { return lambda * phi; }
    virtual ~Polynomial() {}
};

struct Poly_Step : public Potential
{
    Poly_Step( double _m, double _c, double _d, double _phi_step) :
    m{_m}, c{_c}, d{_d}, phi_step{_phi_step} {}
    
    double m;
    double c;
    double d;
    double phi_step;
    
    double V(double phi)   { return 0.5 * m * m * phi * phi * (1 + c * tanh((phi - phi_step) / d)); }
    double dV(double phi)
    {
        return c * m * m * phi * phi * pow(cosh((phi-phi_step) / d), -2) / (2 * d) + m * m * phi * (1 + c * tanh((phi - phi_step) / d));
    }
    double ddV(double phi) { return 0.5 * m * m * ((-2 * c * phi * phi * tanh((phi-phi_step) / d) / pow(d * cosh((phi-phi_step) / d), 2)) + 2 * (c * tanh((phi-phi_step) / d) + 1) + 4 * c * phi / (d * cosh((phi-phi_step) / d)));}
    virtual ~Poly_Step() {}
};

struct Starobinsky : public Potential
{
    Starobinsky( double _m) : m{_m}  {}

    double m;

    double V(double phi)   { return 0.75 * m * m * (1 - exp(-sqrt(2.0/3) * phi)) * (1 - exp(-sqrt(2.0/3) * phi)); }
    double dV(double phi)  { return 0.75 * m * m * sqrt(8.0/3) * exp(-sqrt(2.0/3) * phi) * (1 - exp(-sqrt(2.0/3) * phi)); }
    double ddV(double phi) { return 0.75 * (4.0/3) * m * m * exp(-2*sqrt(2.0/3) * phi) * (2 - exp(sqrt(2.0/3) * phi)); }
    virtual ~Starobinsky() {}
};

struct AxionMonodromy : public Potential
{
    AxionMonodromy(double _m, double _p, double _phi0, double _f0, double _b, double _pf, double _gamma0) : m{_m}, p{_p}, phi0{_phi0}, f0{_f0}, b{_b}, pf{_pf}, gamma0{_gamma0} {}
    
    double m, p, phi0, f0, b, pf, gamma0;
    
    double V(double phi)   { return std::pow(m,4 - p)*(std::pow(std::abs(phi),p) + (b*f0*p*std::pow(phi0,-1 + p)*cos(gamma0 + (std::pow(std::abs(phi)/phi0,1 + pf)*phi0)/f0))/(1 + pf)); }
    double dV(double phi)  { return std::pow(m,4 - p)*p*(std::pow(std::abs(phi),-1 + p) - b*std::pow(std::abs(phi)/phi0,pf)*std::pow(phi0,-1 + p)*sin(gamma0 + (std::abs(phi)*std::pow(std::abs(phi)/phi0,pf))/f0)); }
    double ddV(double phi) { return -(std::pow(m,4 - p)*p*std::pow(std::abs(phi),-2 + p)) + std::pow(m,4 - p)*std::pow(p,2)*std::pow(std::abs(phi),-2 + p) - (b*std::pow(m,4 - p)*p*std::pow(std::abs(phi)/phi0,2*pf)*std::pow(phi0,-1 + p)*cos(gamma0 + (std::abs(phi)*std::pow(std::abs(phi)/phi0,pf))/f0))/f0 - (b*std::pow(m,4 - p)*p*pf*std::pow(std::abs(phi)/phi0,2*pf)*std::pow(phi0,-1 + p)*cos(gamma0 + (std::abs(phi)*std::pow(std::abs(phi)/phi0,pf))/f0))/f0 - b*std::pow(m,4 - p)*p*pf*std::pow(std::abs(phi)/phi0,-1 + pf)*std::pow(phi0,-2 + p)*sin(gamma0 + (std::abs(phi)*std::pow(std::abs(phi)/phi0,pf))/f0); }
    virtual ~AxionMonodromy() {}
};
