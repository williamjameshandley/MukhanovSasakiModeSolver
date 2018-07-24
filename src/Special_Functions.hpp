#pragma once
#include <cmath>
#include <complex>

const std::complex<double> I(0, 1);

std::complex<double> Bessel_J(double v, double x);
std::complex<double> Bessel_Y(double v, double x);

std::complex<double> Hankel1(double v, double x);
std::complex<double> Hankel2(double v, double x);

void Airy(double x, double& Ai, double& Aip, double& Bi, double& Bip);
double Airy_Ai(double x);
double Airy_Bi(double x);
double Airy_Aip(double x);
double Airy_Bip(double x);
double Gamma(double x);
