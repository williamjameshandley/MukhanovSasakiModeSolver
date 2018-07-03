#include <iostream>
#include <Eigen/Dense>
#include <vector>

#include "src/BackgroundSolver.hpp"
#include <boost/numeric/odeint.hpp>


const double sigma = 10.0;
const double R = 28.0;
const double b = 8.0 / 3.0;


void lorenz( const std::vector<double> &x , std::vector<double> &dxdt , double t )
{
    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = -b * x[2] + x[0] * x[1];
}

void write_lorenz( const std::vector<double> &x , const double t )
{
    std::cout << t << " " << x[0] << " " << x[1] << " " << x[2] << std::endl;
}

int main(int argc, char **argv)
{
    std::vector<double> x = { 10.0 , 1.0 , 1.0 }; // initial conditions
    boost::numeric::odeint::integrate( lorenz , x , 0.0 , 25.0 , 0.1 , write_lorenz );
}
