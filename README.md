# MukhanovSasakiModeSolver

## This is a tool for rapid numerical computation of the primordial power spectrum.

##Set up 

* Initialise potential (Choose from Potential.hpp) with parameters and initialise potential pointer. 
* Set initial conditions for background variables: N_star, the number of e-folds before the end of inflation and, N_dagger the number of e-fold between the start of inflation and the reference scale k_*.
* Solve for background variables.
* Choose the vacuum and initial setting time for the Mukhanov-Sasaki equation.
* Initialise ModeSolver with background solutions. 
* Construct primordial power spectrum with required tolerance using a linear interpolator.
* Evaluate the power spectrum from the linear interpolator using ms.PPS(k) or evaluate directly using ms.Find_PPS(k). 



## Potentials 

### Polynomial: 1/2 * m^2 * phi^2 + 1/24 * lambda * phi^4
* Polynomial pot(m)
* Polynomial pot(m, lambda) 

###Stepped Potential: 1/2 * m^2 * phi^2 * (1 + c * tanh((phi - phi_step) / d))
* Poly_Step pot(m, c, d, phi_step)

###Starobinsky: 3/4 * m^2 * (1 - exp(-sqrt(2/3) * phi))^2
* Starobinsky pot(m)
