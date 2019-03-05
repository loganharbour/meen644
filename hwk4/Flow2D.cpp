#include "Flow2D.h"

Flow2D::Flow2D(const unsigned int Nx,
               const unsigned int Ny,
               const double Lx,
               const double Ly,
               const BoundaryCondition u_BC,
               const BoundaryCondition v_BC,
               const double rho,
               const double mu,
               const bool loud,
               const unsigned int max_its,
               const double tol)
  : // Number of pressure CVs
    Nx(Nx),
    Ny(Ny),
    // Maixmum nodal values
    Mx_u(Nx),
    My_u(Ny + 1),
    Mx_v(Nx + 1),
    My_v(Ny),
    Mx_p(Nx + 1),
    My_p(Ny + 1),
    // Sizes
    Lx(Lx),
    Ly(Ly),
    L_ref(Ly),
    dx(Lx / Nx),
    dy(Ly / Ny),
    // Boundary conditions
    u_BC(u_BC),
    v_BC(v_BC),
    // Reference velocity for pressure residual
    u_ref(u_BC.top),
    // Material properties
    rho(rho),
    mu(mu),
    // Enable debug (loud)
    loud(loud),
    // Material properties in matrix form
    a_u(Mx_u + 1, My_u + 1),
    a_v(Mx_v + 1, My_v + 1),
    a_pc(Mx_p + 1, My_p + 1),
    // Solver properties
    max_its(max_its),
    tol(tol),
    w_uv(2.0),
    alpha_p(0.7),
    // Initialize coefficient matrices
    u(Mx_u + 1, My_u + 1),
    v(Mx_v + 1, My_v + 1),
    p(Mx_p + 1, My_p + 1),
    pc(Mx_p + 1, My_p + 1),
    // Initialize sweeping matrices and vectors
    Ax_u(Mx_u - 1),
    Ay_u(My_u - 1),
    Ax_v(Mx_v - 1),
    Ay_v(My_v - 1),
    Ax_pc(Mx_p - 1),
    Ay_pc(My_p - 1),
    bx_u(Mx_u - 1),
    by_u(My_u - 1),
    bx_v(Mx_v - 1),
    by_v(My_v - 1),
    bx_pc(Mx_p - 1),
    by_pc(My_p - 1)
{
}

void
Flow2D::run()
{
  fillInitialValues();

  for (unsigned int l = 0; l < max_its; ++l)
  {
    std::cout << "Iteration " << std::setw(3) << std::left << l << ": ";
    if (loud)
      std::cout << std::endl;

    // Solve for v, u, and then p
    solve(Equations::u);
    solve(Equations::v);
    solve(Equations::pc);

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

  // Correct pressure values on boundary
  pBoundaryCorrect();
}

void
Flow2D::fillInitialValues()
{
  // Set solution values to zero
  u = 0;
  v = 0;
  p = 0;

  // Apply u boundary conditions
  u.setColumn(0, u_BC.left);
  u.setColumn(Mx_u, u_BC.right);
  u.setRow(0, u_BC.bottom);
  u.setRow(My_u, u_BC.top);

  // Apply v boundary conditions
  v.setColumn(0, v_BC.left);
  v.setColumn(Mx_v, v_BC.right);
  v.setRow(0, v_BC.bottom);
  v.setRow(My_v, v_BC.top);
}
