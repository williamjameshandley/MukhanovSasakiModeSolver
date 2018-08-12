#pragma once
#include <cmath>
#include <complex>
#include "cephes.hpp"

const std::complex<double> I(0, 1);

std::complex<double> Bessel_J(double v, double x);
std::complex<double> Bessel_Y(double v, double x);

std::complex<double> Hankel1(double v, double x);
std::complex<double> Hankel2(double v, double x);

void Airy(double x, double& Ai, double& Aip, double& Bi, double& Bip);
double Gamma(double x);
