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
    // Other boundary conditions
    q_top_bot(input.q_top_bot),
    // Mass inflow
    m_in(input.u_bc.left * input.rho * Ly),
    // Material properties
    cp(input.cp),
    k(input.k),
    mu(input.mu),
    rho(input.rho),
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
    T(Variables::T, Nx + 2, Ny + 2, 1),
    // Initialize aux variables
    p(Variables::p, Nx + 2, Ny + 2)
{
  // Add into variable map for access outside of class
  variables.emplace(Variables::u, u);
  variables.emplace(Variables::v, v);
  variables.emplace(Variables::pc, pc);
  variables.emplace(Variables::p, p);
  variables.emplace(Variables::T, T);

  // T initial condition
  T = input.T_left;
}

void
Problem::run()
{
  // Store start time
  start = clock();

  for (unsigned int l = 0; l < max_its; ++l)
  {
    // Solve and correct main variables
    if (!main_converged)
    {
      solveMain();
      correctMain();
      computeMainResiduals();
    }
    // Solve and correct aux variables
    else
    {
      solveAux();
      correctAux();
      computeAuxResiduals();
    }

    // Exit if everything is converged
    if (converged)
      return;
  }

  // Oops. Didn't converge
  cout << "Did not converge after " << main_iterations << " main iterations and " << aux_iterations
       << "aux iterations!" << endl;
}

} // namespace Flow2D
