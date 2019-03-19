#include "Problem.h"

namespace Flow2D
{

Problem::Problem(const unsigned int Nx, const unsigned int Ny, const InputArguments & input)
  : // Number of pressure CVs
    Nx(Nx),
    Ny(Ny),
    // Domain sizes
    Lx(input.Lx),
    Ly(input.Ly),
    dx(Lx / Nx),
    dy(Ly / Ny),
    // Residual references
    L_ref(input.L_ref),
    u_ref(input.u_ref),
    // Material properties
    rho(input.rho),
    mu(input.mu),
    // Enable debug
    debug(input.debug),
    // Solver properties
    max_its(input.max_its),
    tol(input.tol),
    alpha_p(input.alpha_p),
    // Initialize variables for u, v, pc (solved variables)
    u(Variables::u, Nx + 1, Ny + 2, input.alpha_uv, input.u_bc),
    v(Variables::v, Nx + 2, Ny + 1, input.alpha_uv, input.v_bc),
    pc(Variables::pc, Nx + 2, Ny + 2, 1),
    // Initialize aux variables
    p(Variables::p, Nx + 2, Ny + 2)
{
  // Add into variable map for access outside of class
  variables.emplace(Variables::u, u);
  variables.emplace(Variables::v, v);
  variables.emplace(Variables::pc, pc);
  variables.emplace(Variables::p, p);
}

void
Problem::run()
{
  // Store start time
  start = clock();

  for (unsigned int l = 0; l < max_its; ++l)
  {
    ++iterations;
    if (debug)
      cout << "Iteration " << l << endl << endl;

    // Solve for all variables
    solve();

    // Apply corrections
    correct();

    // Compute residuals and exit if converged
    computeResiduals();
    if (converged)
      return;
  }

  // Oops. Didn't converge
  cout << "Did not converge after " << max_its << " iterations!" << endl;
}

}
