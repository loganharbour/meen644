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
Flow2D::solve()
{
  fillInitialValues();

  for (unsigned int l = 0; l < 2; ++l)
  {
    // Reset pressure correction
    pc = 0;

    // Solve for v, u, and then p
    solve(Equations::u);
    solve(Equations::v);
    solve(Equations::pc);

    // Apply corrections
    correct();
    //
    // for (unsigned int j = 0; j <= My_p; ++j)
    // {
    //   p(0, j) = p(1, j);
    //   p(Mx_p, j) = p(Mx_p - 1, j);
    // }
    // for (unsigned int i = 0; i <= Mx_p; ++i)
    // {
    //   p(i, 0) = p(i, 1);
    //   p(i, My_p) = p(i, My_p - 1);
    // }
    std::cout << pResidual() << std::endl;
  }

  u.save("u.csv");
  v.save("v.csv");
  p.save("p.csv");
  pc.save("pc.csv");
}

void
Flow2D::solve(const Equations eq)
{
  unsigned int Mx, My;
  // Fill coefficients for this equation and grab correct max vals
  switch (eq)
  {
    case Equations::u:
      Mx = Mx_u;
      My = My_u;
      uCoefficients();
      break;
    case Equations::v:
      Mx = Mx_v;
      My = My_v;
      vCoefficients();
      break;
    case Equations::pc:
      Mx = Mx_p;
      My = My_p;
      pcCoefficients();
      break;
  }

  for (int i = 1; i < Mx; ++i)
    solveColumn(i, eq);
  for (int j = 1; j < My; ++j)
    solveRow(j, eq);
  // // Solve from bottom to top
  // for (int j = 1; j < My; ++j)
  //   solveRow(j, eq);
  // // // Solve from left to right
  // for (int i = 1; i < Mx; ++i)
  //   solveColumn(i, eq);
  // // Solve from top to bottom
  // for (int j = My - 1; j >= 1; --j)
  //   solveRow(j, eq);
  // // Solve from right to left
  // for (int i = Mx - 1; i >= 1; --i)
  //   solveColumn(i, eq);
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

void
Flow2D::uCoefficients()
{
  Coefficients D, F;
  double W, dy_pn, dy_ps, b;

  for (unsigned int i = 1; i < Mx_u; ++i)
    for (unsigned int j = 1; j < My_u; ++j)
    {
      std::cout << "u coefficients " << i << ", " << j;
      // Width of the cell
      W = (i == 1 || i == Mx_u - 1 ? 3 * dx / 2 : dx);
      // North/south distances to pressure nodes
      dy_pn = (j == My_u - 1 ? dy / 2 : dy);
      dy_ps = (j == 1 ? dy / 2 : dy);

      // Diffusion coefficients
      D.n = mu * W / dy_pn;
      D.e = mu * dy / dx;
      D.s = mu * W / dy_ps;
      D.w = mu * dy / dx;

      // East and west flows
      F.e = (i == Mx_u - 1 ? rho * dy * u(Mx_u, j) : rho * dy * (u(i - 1, j) + u(i, j)) / 2);
      F.w = (i == 1 ? rho * dy * u(0, j) : rho * dy * (u(i - 1, j) + u(i, j)) / 2);
      // North and south flows
      if (i == 1) // Left boundary
      {
        F.n = rho * W * (v(0, j) + 3 * v(1, j) + 2 * v(2, j)) / 6;
        F.s = rho * W * (v(0, j - 1) + 3 * v(1, j - 1) + 2 * v(2, j - 1)) / 6;
      }
      else if (i == Mx_u - 1) // Right boundary
      {
        F.n = rho * W * (2 * v(i, j) + 3 * v(i + 1, j) + v(i + 2, j)) / 6;
        F.s = rho * W * (2 * v(i, j - 1) + 3 * v(i + 1, j - 1) + v(i + 2, j - 1)) / 6;
      }
      else // Interior (not left or right boundary)
      {
        F.n = rho * W * (v(i, j) + v(i + 1, j)) / 2;
        F.s = rho * W * (v(i, j - 1) + v(i + 1, j - 1)) / 2;
      }

      // Pressure RHS
      b = dy * (p(i, j) - p(i + 1, j));

      // Compute Perchlet number and store into a_u(i, j)
      velocityCoefficients(a_u(i, j), D, F, b);
    }
}

void
Flow2D::velocityCoefficients(Coefficients & a,
                                 const Coefficients & D,
                                 const Coefficients & F,
                                 const double b)
{
  // Perchlet number
  Coefficients P;
  P.n = F.n / D.n;
  P.e = F.e / D.e;
  P.s = F.s / D.s;
  P.w = F.w / D.w;

  // Fill coefficients
  a.n = D.n * std::fmax(0, std::pow(1 - 0.1 * std::fabs(P.n), 5)) + std::fmax(-F.n, 0);
  a.e = D.e * std::fmax(0, std::pow(1 - 0.1 * std::fabs(P.e), 5)) + std::fmax(-F.e, 0);
  a.s = D.s * std::fmax(0, std::pow(1 - 0.1 * std::fabs(P.s), 5)) + std::fmax(F.s, 0);
  a.w = D.w * std::fmax(0, std::pow(1 - 0.1 * std::fabs(P.w), 5)) + std::fmax(F.w, 0);
  a.p = a.n + a.e + a.s + a.w;
  a.b = b;

  std::cout << ":  n = " << a.n;
  std::cout << ", e = " << a.e;
  std::cout << ", s = " << a.s;
  std::cout << ", w = " << a.w;
  std::cout << ", b = " << a.b << std::endl;
}

void
Flow2D::vCoefficients()
{
  Coefficients D, F;
  double H, dx_pe, dx_pw, b;

  for (unsigned int i = 1; i < Mx_v; ++i)
    for (unsigned int j = 1; j < My_v; ++j)
    {
      std::cout << "v coefficients " << i << ", " << j;
      // Height of the cell
      H = (j == 1 || j == My_v - 1 ? 3 * dy / 2 : dy);
      // East/west distances to pressure nodes
      dx_pe = (i == Mx_v - 1 ? dx / 2 : dx);
      dx_pw = (i == 1 ? dx / 2 : dx);

      // Diffusion coefficient
      D.n = mu * dx / dy;
      D.e = mu * H / dx_pe;
      D.s = mu * dx / dy;
      D.w = mu * H / dx_pw;

      // North and east flows
      F.n = (j == My_v - 1 ? rho * dx * v(i, My_v) : rho * dx * (v(i, j + 1) + v(i, j)) / 2);
      F.s = (j == 1 ? rho * dx * v(i, 0) : rho * dx * (v(i, j - 1) + v(i, j)) / 2);
      // East and west flows
      if (j == 1) // Bottom boundary
      {
        F.e = rho * H * (u(i, 0) + 3 * u(i, 1) + 2 * u(i, 2)) / 6;
        F.w = rho * H * (u(i - 1, 0) + 3 * u(i - 1, 1) + 2 * u(i - 1, 2)) / 6;
      }
      else if (j == My_v - 1) // Top boundary
      {
        F.e = rho * H * (u(i, My_u) + 3 * u(i, My_u - 1) + 2 * u(i, My_u - 2)) / 6;
        F.w = rho * H * (u(i - 1, My_u) + 3 * u(i - 1, My_u - 1) + 2 * u(i - 1, My_u - 2)) / 6;
      }
      else // Interior (not top or bottom boundary)
      {
        F.e = rho * H * (u(i, j) + u(i, j + 1)) / 2;
        F.w = rho * H * (u(i - 1, j) + u(i - 1, j + 1)) / 2;
      }

      // Pressure RHS
      b = dx * (p(i, j) - p(i, j + 1));

      // Compute Perchlet number and store into a_v(i, j)
      velocityCoefficients(a_v(i, j), D, F, b);
    }
}

void
Flow2D::solveRow(const unsigned int j, const Equations eq)
{
  // Grab the appropriate matrix/vector/coefficients/previous solution
  auto & A = (eq == Equations::u ? Ax_u : (eq == Equations::v ? Ax_v : Ax_pc));
  auto & b = (eq == Equations::u ? bx_u : (eq == Equations::v ? bx_v : bx_pc));
  auto M = (eq == Equations::u ? Mx_u : (eq == Equations::v ? Mx_v : Mx_p));
  Matrix<double> & phi = (eq == Equations::u ? u : (eq == Equations::v ? v : pc));
  auto & as = (eq == Equations::u ? a_u : (eq == Equations::v ? a_v : a_pc));
  double w = (eq == Equations::pc ? 1 : w_uv);

  // Fill for each cell
  for (unsigned int i = 1; i < M; ++i)
  {
    Coefficients & a = as(i, j);
    b[i - 1] = a.b + a.s * phi(i, j - 1) + a.n * phi(i, j + 1) + a.p * phi(i, j) * (w - 1);
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
  A.save("row_" + std::to_string(j) + "_A.csv");
  saveCSV(b, "row_" + std::to_string(j) + "_b.csv");
  A.solveTDMA(b);
  saveCSV(b, "row_" + std::to_string(j) + "_sol.csv");
  for (unsigned int i = 1; i < M; ++i)
    phi(i, j) = b[i - 1];
}

void
Flow2D::solveColumn(const unsigned int i, const Equations eq)
{
  // Grab the appropriate matrix/vector/coefficients/previous solution
  auto & A = (eq == Equations::u ? Ay_u : (eq == Equations::v ? Ay_v : Ay_pc));
  auto & b = (eq == Equations::u ? by_u : (eq == Equations::v ? by_v : by_pc));
  auto M = (eq == Equations::u ? My_u : (eq == Equations::v ? My_v : My_p));
  auto & phi = (eq == Equations::u ? u : (eq == Equations::v ? v : pc));
  auto & as = (eq == Equations::u ? a_u : (eq == Equations::v ? a_v : a_pc));
  double w = (eq == Equations::pc ? 1 : w_uv);

  // Fill for each cell
  for (unsigned int j = 1; j < M; ++j)
  {
    Coefficients & a = as(i, j);
    b[j - 1] = a.b + a.w * phi(i - 1, j) + a.e * phi(i + 1, j) + a.p * phi(i, j) * (w - 1);
    if (j == 1)
    {
      A.setTopRow(a.p * w, -a.n);
      if (eq != Equations::pc)
        b[j - 1] += a.s * phi(i, j - 1);
    }
    else if (j == M - 1)
    {
      A.setBottomRow(-a.s, a.p * w);
      if (eq != Equations::pc)
        b[j - 1] += a.n * phi(i, j + 1);
    }
    else
      A.setMiddleRow(j - 1, -a.s, a.p * w, -a.n);
  }

  // And solve
  A.save("col_" + std::to_string(i) + "_A.csv");
  saveCSV(b, "col_" + std::to_string(i) + "_b.csv");
  A.solveTDMA(b);
  saveCSV(b, "col_" + std::to_string(i) + "_sol.csv");
  for (unsigned int j = 1; j < M; ++j)
    phi(i, j) = b[j - 1];
}

void
Flow2D::pcCoefficients()
{
  for (unsigned int i = 1; i < Mx_p; ++i)
    for (unsigned int j = 1; j < My_p; ++j)
    {
      Coefficients & a = a_pc(i, j);

      // Each of the above is only valid off boundary
      if (i != 1)
        a.w = rho * dy * dy / a_u(i - 1, j).p;
      if (i != Mx_p - 1)
        a.e = rho * dy * dy / a_u(i, j).p;
      if (j != 1)
        a.s = rho * dx * dx / a_v(i, j - 1).p;
      if (j != My_p - 1)
        a.n = rho * dx * dx / a_v(i, j).p;
      a.p = a.n + a.e + a.s + a.w;
      a.b = rho * (dy * (u(i - 1, j) - u(i, j)) + dx * (v(i, j - 1) + v(i, j)));
      std::cout << "(" << i << ", " << j << ")";
      std::cout << " an = " << a.n;
      std::cout << ", ae = " << a.e;
      std::cout << ", as = " << a.s;
      std::cout << ", aw = " << a.w;
      std::cout << ", ab = " << a.b;
      std::cout << ", avp = " << a_v(i, j).p << std::endl;
    }
}

void
Flow2D::correct()
{
  // Correct u-velocities
  for (unsigned int i = 1; i < Mx_u; ++i)
    for (unsigned int j = 1; j < My_u; ++j)
      u(i, j) += dy * (pc(i, j) - pc(i + 1, j)) / a_u(i, j).p;

  // Correct v-velocities
  for (unsigned int i = 1; i < Mx_v; ++i)
    for (unsigned int j = 1; j < My_v; ++j)
      v(i, j) += dx * (pc(i, j) - pc(i, j + 1)) / a_v(i, j).p;

  // Correct pressures
  for (unsigned int i = 1; i < Mx_p; ++i)
    for (unsigned int j = 1; j < My_p; ++j)
      p(i, j) += alpha_p * pc(i, j);
}

double Flow2D::pResidual() {
  double R = 0;
  for (unsigned int i = 1; i < Mx_p; ++i)
    for (unsigned int j = 1; i < My_p; ++j)
      P += std::fabs(rho * (dy * (u(i - 1, j) - u(i, j)) + dx * (v(i, j - 1) - v(i, j))) / (Re * mu);
  return P;
}
