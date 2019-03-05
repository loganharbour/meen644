#include "Flow2D.h"

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
      a.b = rho * (dy * (u(i - 1, j) - u(i, j)) + dx * (v(i, j - 1) - v(i, j)));
    }

  if (loud)
  {
    std::cout << "pc coefficients: " << std::endl;
    a_pc.print();
    std::cout << std::endl;
  }
}

void
Flow2D::uCoefficients()
{
  Coefficients D, F;
  double W, dy_pn, dy_ps, b;

  for (unsigned int i = 1; i < Mx_u; ++i)
    for (unsigned int j = 1; j < My_u; ++j)
    {
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
      F.e = (i == Mx_u - 1 ? rho * dy * u(Mx_u, j) : rho * dy * (u(i + 1, j) + u(i, j)) / 2);
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

  if (loud)
  {
    std::cout << "u coefficients: " << std::endl;
    a_u.print();
    std::cout << std::endl;
  }
}

void
Flow2D::vCoefficients()
{
  Coefficients D, F;
  double H, dx_pe, dx_pw, b;

  for (unsigned int i = 1; i < Mx_v; ++i)
    for (unsigned int j = 1; j < My_v; ++j)
    {
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

  if (loud)
  {
    std::cout << "v coefficients: " << std::endl;
    a_v.print();
    std::cout << std::endl;
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
}
