#pragma once
#include <cmath>
#include <complex>
#include "cephes.hpp"

const std::complex<double> I(0, 1);

std::complex<double> Bessel_J(int v, double x);
std::complex<double> Bessel_Y(int v, double x);

std::complex<double> Bessel_I(int v, double x);
std::complex<double> Bessel_K(int v, double x);

std::complex<double> Hankel1(int v, double x);
std::complex<double> Hankel2(int v, double x);

void Airy(double x, double& Ai, double& Aip, double& Bi, double& Bip);

