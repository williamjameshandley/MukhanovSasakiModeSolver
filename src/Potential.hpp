#pragma once
#include <iostream>

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

