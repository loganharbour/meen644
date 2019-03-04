#include "Flow2D.h"

#include <cmath>

Flow2D::Flow2D(unsigned int Nx,
               unsigned int Ny,
               double Lx,
               double Ly,
               BoundaryCondition u_BC,
               BoundaryCondition v_BC,
               double rho,
               double mu,
               unsigned int max_its)
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
    dx(Lx / Nx),
    dy(Ly / Ny),
    // Boundary conditions
    u_BC(u_BC),
    v_BC(v_BC),
    // Material properties
    rho(rho),
    mu(mu),
    // Material properties in matrix form
    a_u(Mx_u + 1, My_u + 1),
    a_v(Mx_v + 1, My_v + 1),
    a_pc(Mx_p + 1, My_p + 1),
    // Solver properties
    max_its(max_its),
    w_uv(2.0),
    alpha_p(0.7),
    // Initialize coefficient matrices
    u(Mx_u + 1, My_u + 1),
    v(Mx_v + 1, My_v + 1),
    p(Mx_p + 1, My_p + 1),
    pc(Mx_p + 1, My_p + 1),
    // Initialize sweeping matrices and vectors
    Ax_u(My_u - 1),
    Ay_u(Mx_u - 1),
    Ax_v(My_v - 1),
    Ay_v(Mx_v - 1),
    Ax_pc(My_p - 1),
    Ay_pc(Mx_p - 1),
    bx_u(My_u - 1),
    by_u(Mx_u - 1),
    bx_v(My_v - 1),
    by_v(Mx_v - 1),
    bx_pc(My_p - 1),
    by_pc(Mx_p - 1)
{
}

void
Flow2D::solve()
{
  // for (unsigned int l = 0; l < max_its; ++l)
  // {
  // Solve for v, u, and then p
  solve(Equations::u);
  solve(Equations::v);
  solve(Equations::pc);

  // Apply corrections
  // correct();
  // }

  u.save("u.csv");
  v.save("v.csv");
  p.save("p.csv");
  pc.save("pc.csv");
}

void
Flow2D::solve(Equations eq)
{
  unsigned int Mx, My;
  // Fill coefficients for this equation and grab correct max vals
  switch (eq)
  {
    case Equations::u:
      Mx = Mx_u;
      My = My_u;
      filluCoefficients();
      break;
    case Equations::v:
      Mx = Mx_v;
      My = My_v;
      fillvCoefficients();
      break;
    case Equations::pc:
      Mx = Mx_p;
      My = My_p;
      fillpcCoefficients();
      break;
  }

  // Solve from bottom to top
  for (int j = 1; j < My; ++j)
    solveRow(j, eq);
  // Solve from left to right
  for (int i = 1; i < Mx; ++i)
    solveColumn(i, eq);
  // Solve from top to bottom
  for (int j = My - 1; j >= 1; --j)
    solveRow(j, eq);
  // Solve from right to left
  for (int i = Mx - 1; i >= 1; --i)
    solveColumn(i, eq);
}

void
Flow2D::fillBCs()
{
  u.setRow(0, u_BC.bottom);
  u.setRow(My_u, u_BC.top);
  u.setColumn(0, u_BC.left);
  u.setColumn(Mx_u, u_BC.right);
  v.setRow(0, v_BC.bottom);
  v.setRow(My_v, v_BC.top);
  v.setColumn(0, v_BC.left);
  v.setColumn(Mx_v, v_BC.right);
}

void
Flow2D::filluCoefficients()
{
  Coefficients D, F, P;

  for (unsigned int i = 1; i < Mx_u; ++i)
    for (unsigned int j = 1; j < My_u; ++j)
    {
      // Diffusion coefficient for left and right
      D.e = dy * mu / dx;
      D.w = dy * mu / dx;
      // Diffusion coefficient for top and bottom for internal cells
      if (i > 1 && i < Mx_u - 1)
      {
        D.n = 3 * dx * mu / (2 * dy);
        D.s = 3 * dx * mu / (2 * dy);
      }
      // Diffusion coefficient for top and bottom for left and right cells
      else
      {
        D.n = dx * mu / dy;
        D.s = dx * mu / dy;
      }
      if (j == 1)
        D.s *= 2;
      if (j == My_u - 1)
        D.n *= 2;

      // West flow rates
      if (i == 1)
        F.w = rho * dy * u(0, j);
      else
        F.w = rho * dy * (u(i - 1, j) + u(i, j)) / 2;
      // East flow rates
      if (i == Mx_u - 1)
        F.w = rho * dy * u(Mx_u, j);
      else
        F.w = rho * dy * (u(i - 1, j) + u(i, j)) / 2;
      // North and south flow rates on left boundary
      if (i == 1)
      {
        F.n = rho * dx * (v(0, j) + 3 * v(1, j) + 2 * v(2, j)) / 4;
        F.s = rho * dx * (v(0, j - 1) + 3 * v(1, j - 1) + 2 * v(2, j - 1)) / 4;
      }
      // North and south flow rates on right boundary
      else if (i == Mx_u - 1)
      {
        F.n = rho * dx * (2 * v(i, j) + 3 * v(i + 1, j) + v(i + 2, j)) / 4;
        F.s = rho * dx * (2 * v(i, j - 1) + 3 * v(i + 1, j - 1) + v(i + 2, j - 1)) / 4;
      }
      // North and south flow rates on the remainder
      else
      {
        F.n = rho * dx * (v(i, j) + v(i + 1, j)) / 2;
        F.s = rho * dx * (v(i, j - 1) + v(i + 1, j - 1)) / 2;
      }

      // Perchlet number
      P.n = F.n / D.n;
      P.e = F.e / D.e;
      P.s = F.s / D.s;
      P.w = F.w / D.w;
      std::cout << "(" << i << "," << j << "): n = " << D.n << ", e " << D.e << ", s " << D.s << ", w " << D.w << std::endl;

      // Fill coefficients
      Coefficients & a = a_u(i, j);
      a.n = D.n * std::fmax(0, std::pow(1 - 0.1 * std::fabs(P.n), 5)) + std::fmax(-F.n, 0);
      a.e = D.e * std::fmax(0, std::pow(1 - 0.1 * std::fabs(P.e), 5)) + std::fmax(-F.e, 0);
      a.s = D.s * std::fmax(0, std::pow(1 - 0.1 * std::fabs(P.s), 5)) + std::fmax(F.s, 0);
      a.w = D.w * std::fmax(0, std::pow(1 - 0.1 * std::fabs(P.w), 5)) + std::fmax(F.w, 0);
      a.p = a.n + a.e + a.s + a.w;
      // std::cout << "(" << i << "," << j << "): n = " << a.n << ", e " << a.e << ", s " << a.s << ", w " << a.w << ", p " << a.p << std::endl;
      a.b = dy * (p(i, j) - p(i + 1, j));
    }
}

void
Flow2D::fillvCoefficients()
{
  Coefficients D, F, P;

  for (unsigned int i = 1; i < Mx_v; ++i)
    for (unsigned int j = 1; j < My_v; ++j)
    {
      // Diffusion coefficient for top and bottom
      D.n = dx * mu / dy;
      D.s = dx * mu / dy;
      // Diffusion coefficient for left and right for internal cells
      if (j > 1 && j < Mx_v - 1)
      {
        D.e = 3 * dy * mu / (2 * dx);
        D.w = 3 * dy * mu / (2 * dx);
      }
      // Diffusion coefficient for left and right for bottom and top cells
      else
      {
        D.e = dy * mu / dx;
        D.w = dy * mu / dx;
      }
      // Why
      if (i == 1)
        D.w *= 2;
      if (i == Mx_v - 1)
        D.e *= 2;

      // North flow rates
      if (j == My_v - 1)
        F.n = rho * dx * v(i, My_v);
      else
        F.n = rho * dx * (v(i, j + 1) + v(i, j)) / 2;
      // South flow rates
      if (j == 1)
        F.s = rho * dx * v(i, 0);
      else
        F.s = rho * dx * (v(i, j - 1) + v(i, j)) / 2;
      // East and west flow rates on bottom boundary
      if (j == 1)
      {
        F.e = rho * dy * (u(i, 0) + 2 * u(i, 1) + 3 * u(i, 2)) / 4;
        F.w = rho * dy * (u(i - 1, 0) + 2 * u(i - 1, 1) + 3 * u(i - 1, 2)) / 4;
      }
      // East and west flow rates on top boundary
      else if (j == My_v - 1)
      {
        F.e = rho * dy * (u(i, My_u) + 2 * u(i, My_u - 1) + 3 * u(i, My_u - 2)) / 4;
        F.w = rho * dy * (u(i - 1, My_u) + 2 * u(i - 1, My_u - 1) + 3 * u(i - 1, My_u - 2)) / 4;
      }
      // East and west flow rates on the remainder
      else
      {
        F.e = rho * dy * (u(i, j) + u(i, j + 1)) / 2;
        F.w = rho * dy * (u(i - 1, j) + u(i - 1, j + 1)) / 2;
      }

      // Perchlet number
      P.n = F.n / D.n;
      P.e = F.e / D.e;
      P.s = F.s / D.s;
      P.w = F.w / D.w;

      // Fill coefficients
      Coefficients & a = a_v(i, j);
      a.n = D.n * std::fmax(0, std::pow(1 - 0.1 * std::fabs(P.n), 5)) + std::fmax(-F.n, 0);
      a.e = D.e * std::fmax(0, std::pow(1 - 0.1 * std::fabs(P.e), 5)) + std::fmax(-F.e, 0);
      a.s = D.s * std::fmax(0, std::pow(1 - 0.1 * std::fabs(P.s), 5)) + std::fmax(F.s, 0);
      a.w = D.w * std::fmax(0, std::pow(1 - 0.1 * std::fabs(P.w), 5)) + std::fmax(F.w, 0);
      a.p = a.n + a.e + a.s + a.w;
      a.b = dx * (p(i, j) - p(i + 1, j));
    }
}

void
Flow2D::solveRow(unsigned int j, Equations eq)
{
  // Grab the appropriate matrix/vector/coefficients/previous solution
  auto & A = (eq == Equations::u ? Ax_u : (eq == Equations::v ? Ax_u : Ax_pc));
  auto & b = (eq == Equations::u ? bx_u : (eq == Equations::v ? bx_u : bx_pc));
  auto M = (eq == Equations::u ? Mx_u : (eq == Equations::v ? Mx_u : Mx_p));
  auto & phi = (eq == Equations::u ? u : (eq == Equations::v ? v : pc));
  auto & as = (eq == Equations::u ? a_u : (eq == Equations::v ? a_v : a_pc));
  double w = (eq == Equations::pc ? 1 : w_uv);

  // Fill for each cell
  for (unsigned int i = 1; i < M; ++i)
  {
    Coefficients & a = as(i, j);
    b[i - 1] = a.b + a.s * phi(i, j - 1) + a.n * phi(i, j + 1);
    if (w != 1)
      b[i - 1] += a.p * phi(i, j) * (w - 1);
    if (i == 1)
    {
      A.setTopRow(a.p * w, -a.e);
      if (eq != Equations::pc)
        b[i - 1] += a.w * phi(i - 1, j);
    }
    else if (i == M - 1)
    {
      A.setBottomRow(-a.w, a.p * w);
      if (eq != Equations::pc)
        b[i - 1] += a.e * phi(i + 1, j);
    }
    else
      A.setMiddleRow(i - 1, -a.w, a.p * w, -a.e);
  }

  // And solve
  A.solveTDMA(b);
  for (unsigned int i = 1; i < M; ++i)
    phi(i, j) = b[i - 1];
}

void
Flow2D::solveColumn(unsigned int i, Equations eq)
{
  // Grab the appropriate matrix/vector/coefficients/previous solution
  auto & A = (eq == Equations::u ? Ay_u : (eq == Equations::v ? Ay_u : Ay_pc));
  auto & b = (eq == Equations::u ? by_u : (eq == Equations::v ? by_u : by_pc));
  auto M = (eq == Equations::u ? My_u : (eq == Equations::v ? My_u : My_p));
  auto & phi = (eq == Equations::u ? u : (eq == Equations::v ? v : pc));
  auto & as = (eq == Equations::u ? a_u : (eq == Equations::v ? a_v : a_pc));
  double w = (eq == Equations::pc ? 1 : w_uv);

  // Fill for each cell
  for (unsigned int j = 1; j < M; ++j)
  {
    Coefficients & a = as(i, j);
    b[i - 1] = a.b + a.w * phi(i - 1, j) + a.e * phi(i + 1, j);
    if (w != 1)
      b[i - 1] += a.p * phi(i, j) * (w - 1);
    if (j == 1)
    {
      A.setTopRow(a.p * w, -a.n);
      if (eq != Equations::pc)
        b[i - 1] += a.s * phi(i, j - 1);
    }
    else if (j == M)
    {
      A.setBottomRow(-a.s, a.p * w);
      if (eq != Equations::pc)
        b[i - 1] += a.n * phi(i, j + 1);
    }
    else
      A.setMiddleRow(i - 1, -a.s * w, a.p, -a.n);
  }

  // And solve
  A.solveTDMA(b);
  for (unsigned int j = 1; j < M; ++j)
    u(i, j) = b[i - 1];
}

void
Flow2D::fillpcCoefficients()
{
  for (unsigned int i = 1; i < Mx_p; ++i)
    for (unsigned int j = 1; j < My_p; ++j)
    {
      Coefficients & a = a_pc(i, j);

      // Each of the above is only valid off boundary
      if (i != 1)
        a.w = rho * dy * dy / a_v(i - 1, j).p;
      if (i != Mx_p - 1)
        a.e = rho * dy * dy / a_v(i, i).p;
      if (j != 1)
        a.s = rho * dx * dx / a_v(i, j - 1).p;
      if (j != My_p - 1)
        a.n = rho * dx * dx / a_v(i, j).p;
      a.p = a.n + a.e + a.s + a.w;
      a.b = rho * (dy * (u(i - 1, j) - u(i, j)) + dx * (v(i, j - 1) + v(i, j)));
    }
}

void
Flow2D::correct()
{
  // Correct u-velocities
  for (unsigned int i = 1; i < Mx_u; ++i)
    for (unsigned int j = 1; j < My_u; ++j)
      u(i, j) += dy * (pc(i, j) - pc(i + 1, j)) * dy / a_u(i, j).p;

  // Correct v-velocities
  for (unsigned int i = 1; i < Mx_v; ++i)
    for (unsigned int j = 1; j < My_v; ++j)
      v(i, j) += dx * (pc(i, j) - pc(i, j + 1)) * dy / a_v(i, j).p;

  // Correct pressures
  for (unsigned int i = 1; i < Mx_p; ++i)
    for (unsigned int j = 1; j < My_p; ++j)
      p(i, j) += alpha_p * pc(i, j);
}
