#pragma once
#include <iostream>
#include <math.h>

class Poly
{
    public:
        
        double m;
        double lambda;
        
        double V(double phi)
        {
            double V = 0.5 * m * m * phi * phi + (1.0/24) * lambda * phi * phi * phi * phi;
            return V;
        }
        double dV(double phi)
        {
            double dV = m * m * phi + (1.0/6) * lambda * phi * phi * phi;
            return dV;
        }
        double ddV(double phi)
        {
            double ddV = m * m + (1.0/2) * lambda * phi * phi;
            return ddV;
        }
};

class R2
{
public:
    
    double m;
    
    double V(double phi)
    {
        double V = 0.75 * m * m * (1 - exp(-sqrt(2.0/3) * phi)) * (1 - exp(-sqrt(2.0/3) * phi));
        return V;
    }
    double dV(double phi)
    {
        double dV = 0.75 * m * m * sqrt(8.0/3) * exp(-sqrt(2.0/3) * phi) * (1 - exp(-sqrt(2.0/3) * phi));
        return dV;
    }
    double ddV(double phi)
    {
        double ddV = 0.75 * (4.0/3) * m * m * exp(-2*sqrt(2.0/3) * phi) * (2 - exp(sqrt(2.0/3) * phi));
        return ddV;
    }
};

