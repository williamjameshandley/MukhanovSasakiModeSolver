#include <iostream>
#include <fstream>
#include <odepack/dlsodar.hpp>
#include <cmath>
#define USR_FULL_JAC 1

void simple_pendulum_field(double qdot[], const double t, const double q[], double params[]);
void simple_pendulum_jacobian(double dfdq[], const double t, const double q[], double params[]);
void simple_pendulum_cycle(double g[], const double t, const double q[], double params[]);

int main(){
    dlsodar solver(2, 1);

    double q[2] = {M_PI * 99999/100000.0, 0.};
    auto t = 0.0;
    auto tf = 1000.0;
    auto dt = 0.1;
    double alpha = 1.0;
    double params[2] = {alpha, 0};

    std::ofstream fout{"simple_pendulum_trajectory"};
    solver.integrate(
            t, tf, q, 
            simple_pendulum_field, 
            simple_pendulum_cycle, 
            params);
    std::cout << t << " " << q[0] << " " << q[1] << std::endl;

    while(t < tf){    
        solver.integrate(
                t, t+dt, q, 
                simple_pendulum_field, 
                simple_pendulum_cycle, 
                params);
        fout << t << " " << q[0] << " " << q[1] << std::endl;
    }

}

void simple_pendulum_field(double qdot[], const double t, const double q[], double params[])
{
    auto alpha = params[0];

    qdot[0] = q[1];
    qdot[1] = - alpha * sin(q[0]);
}

void simple_pendulum_jacobian(double dfdq[], const double t, const double q[], double params[])
{
    auto alpha = params[0];

    /*Column ordered; ODEPACK is in FORTRAN 77.*/
    dfdq[0] = 0.0;                 dfdq[2] = 1.0;
    dfdq[1] = - alpha * cos(q[0]); dfdq[3] = 0.0;
    std::cout << "Jacobian " << dfdq[1] << std::endl;
}

void simple_pendulum_cycle(double g[], const double t, const double q[], double params[])
{
    auto alpha = params[0];
    auto theta_0 = params[1];

    g[0] = (q[0] - theta_0);
}
