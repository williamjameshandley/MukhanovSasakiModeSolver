MukhanovSasakiModeSolver
========================

This is a tool for rapid numerical computation of the primordial power spectrum.
--------------------------------------------------------------------------------

Set up:

-   Initialise potential (Choose from Potential.hpp) with parameters and
    initialise potential pointer.

-   Set initial conditions for background variables: N_star, the number of
    e-folds before the end of inflation and, N_dagger the number of e-fold
    between the start of inflation and the reference scale k_\*.

-   Solve for background variables

\`\`\`bash

auto sols = solve_equations(Error, potential_ptr, N_star, N_dagger);

\`\`\`

-   Initialise ModeSolver with background solutions

\`\`\`bash

ModeSolver ms(sols);

\`\`\`

-   Choose the vacuum and initial setting time (no. e-folds before end of
    inflation) for the Mukhanov-Sasaki equation.

\`\`\`bash

ms.Initial_Conditions(BD, N_r);

\`\`\`

-   Evaluate the power spectrum directly using

\`\`\`bash

ms.Find_PPS_Scalar(k)

ms.Find_PPS_Tensor(k)

\`\`\`
