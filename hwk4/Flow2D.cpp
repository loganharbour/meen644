#include "Flow2D.h"

Flow2D::Flow2D(const unsigned int Nx,
               const unsigned int Ny,
               const double Lx,
               const double Ly,
               const BoundaryCondition u_bc,
               const BoundaryCondition v_bc,
               const double rho,
               const double mu,
               const bool loud,
               const unsigned int max_its,
               const double tol)
  : // Number of pressure CVs
    Nx(Nx),
    Ny(Ny),
    // Maixmum nodal values
    // Sizes
    Lx(Lx),
    Ly(Ly),
    L_ref(Ly),
    dx(Lx / Nx),
    dy(Ly / Ny),
    // Reference velocity for pressure residual
    u_ref(u_BC.top),
    // Material properties
    rho(rho),
    mu(mu),
    // Enable debug (loud)
    loud(loud),
    // Solver properties
    max_its(max_its),
    tol(tol),
    w_uv(2.0),
    alpha_p(0.7),
    // Initialize systems for u, v, pc
    u(Equations::u, Nx + 1, Ny + 2, u_BC),
    v(Equations::v, Nx + 2, Ny + 1, v_BC),
    pc(Equations::pc, Nx + 2, Ny + 2)
    // Doesn't have a system and is an aux variable
    p(Nx + 2, Ny + 2),
{
}

void
Flow2D::run()
{
  for (unsigned int l = 0; l < max_its; ++l)
  {
    std::cout << "Iteration " << std::setw(3) << std::left << l << ": ";
    if (loud)
      std::cout << std::endl;

    // Solve for u, v, and then p
    solve(u);
    solve(v);
    solve(pc);

    // Apply corrections
    correct();

    // Compute residuals
    computeResiduals();

    // Exit if converged
    if (converged)
      break;
  }

  if (!converged)
    std::cout << "Did not converge!" << std::endl;
}
