#include "Problem.h"

namespace Flow2D
{

Problem::Problem(const unsigned int Nx,
                 const unsigned int Ny,
                 const InputArguments & input)
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
    // Heat flux boundary condition
    q(input.q),
    // Initialize material properties
    k(Nx + 1, Ny + 1),
    mu_u(Nx, Ny + 1),
    mu_v(Nx + 1, Ny),
    rho(input.rho),
    cp(input.cp),
    // Enable debug
    debug(input.debug),
    // Solver properties
    max_main_its(input.max_main_its),
    max_aux_its(input.max_aux_its),
    tol(input.tol),
    alpha_p(input.alpha_p),
    // Initialize variables for u, v, pc (solved variables)
    u(Variables::u, Nx + 1, Ny + 2, dx, dy, input.alpha_uv, input.u_ic),
    v(Variables::v, Nx + 2, Ny + 1, dx, dy, input.alpha_uv, input.v_ic),
    pc(Variables::pc, Nx + 2, Ny + 2, dx, dy, 1),
    T(Variables::T, Nx + 2, Ny + 2, dx, dy, 1, input.T_ic),
    // Initialize aux variables
    p(Variables::p, Nx + 2, Ny + 2, dx, dy)
{
  // Add into variable map for access outside of class
  variables.emplace(Variables::u, u);
  variables.emplace(Variables::v, v);
  variables.emplace(Variables::pc, pc);
  variables.emplace(Variables::p, p);
  variables.emplace(Variables::T, T);

  // Fill non-constant materials
  fillMaterial(k, input.k, T);
  fillMaterial(mu_u, input.mu, u);
  fillMaterial(mu_v, input.mu, v);
}

void
Problem::run()
{
  // Store start time
  start = clock();

  // Solve main variables
  for (unsigned int l = 0; l < max_main_its; ++l)
  {
    solveMain();
    correctMain();
    computeMainResiduals();

    // Break out if we've converged
    if (main_converged)
      break;
  }

  // Ensure main variables converged
  if (!main_converged)
    cout << "Main variables did not converge after " << max_main_its << " iterations!" << endl;

  // Solve aux variables
  for (unsigned int l = 0; l < max_aux_its; ++l)
  {
    solveAux();
    correctAux();
    computeAuxResiduals();

    // Exit if everything is converged
    if (converged)
      return;
  }

  // Oops. Didn't converge
  cout << "Aux variables did not converge after " << max_aux_its << " iterations!" << endl;
}

void
Problem::fillMaterial(Matrix<Coefficients> & m,
                      std::function<double(const std::vector<double> &)> func,
                      const Variable & var)
{
  // First, fill the variable everywhere (p, n, e, s, w)
  for (unsigned int i = 1; i < var.Mx; ++i)
    for (unsigned int j = 1; j < var.My; ++j)
      m(i, j).p = func(var.point(i, j));

  // And now fill with the harmonic mean at the interior edges
  for (unsigned int i = 1; i < var.Mx; ++i)
    for (unsigned int j = 1; j < var.My; ++j)
    {
      if (j != var.My - 1)
        m(i, j).n = 2 * m(i, j).p * m(i, j + 1).p / (m(i, j).p + m(i, j + 1).p);
      else
        m(i, j).n = m(i, j).p;
      if (i != var.Mx - 1)
        m(i, j).e = 2 * m(i, j).p * m(i + 1, j).p / (m(i, j).p + m(i + 1, j).p);
      else
        m(i, j).e = m(i, j).p;
      if (j != 1)
        m(i, j).s = 2 * m(i, j).p * m(i, j - 1).p / (m(i, j).p + m(i, j - 1).p);
      else
        m(i, j).s = m(i, j).p;
      if (i != 1)
        m(i, j).w = 2 * m(i, j).p * m(i - 1, j).p / (m(i, j).p + m(i - 1, j).p);
      else
        m(i, j).w = m(i, j).p;
    }
}
} // namespace Flow2D
