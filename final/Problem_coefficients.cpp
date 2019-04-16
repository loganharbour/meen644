#include "Problem.h"

namespace Flow2D
{

void
Problem::fillCoefficients(const Variable & var)
{
  if (var.name == Variables::pc)
    pcCoefficients();
  else if (var.name == Variables::T)
    TCoefficients();
  else if (var.name == Variables::u)
    uCoefficients();
  else if (var.name == Variables::v)
    vCoefficients();

  if (debug)
  {
    cout << var.string << " coefficients: " << endl;
    var.printCoefficients(var.string, true);
  }
}

void
Problem::pcCoefficients()
{
  for (unsigned int i = 1; i < pc.Mx; ++i)
    for (unsigned int j = 1; j < pc.My; ++j)
    {
      Coefficients & a = pc.a(i, j);

      if (i != 1)
        a.w = rho * dy * dy / u.a(i - 1, j).p;
      if (i != pc.Mx - 1)
        a.e = rho * dy * dy / u.a(i, j).p;
      if (j != 1)
        a.s = rho * dx * dx / v.a(i, j - 1).p;
      if (j != pc.My - 1)
        a.n = rho * dx * dx / v.a(i, j).p;
      a.p = a.n + a.e + a.s + a.w;
      a.b = rho * (dy * (u(i - 1, j) - u(i, j)) + dx * (v(i, j - 1) - v(i, j)));
    }
}

void
Problem::TCoefficients()
{
  for (unsigned int i = 1; i < T.Mx; ++i)
    for (unsigned int j = 1; j < T.My; ++j)
    {
      Coefficients D, F;

      // Diffusion coefficient
      D.n = k(i, j).n * dx / dy * (j == T.My - 1 ? 2.0 : 1.0);
      D.e = k(i, j).e * dy / dx * (i == T.Mx - 1 ? 2.0 : 1.0);
      D.s = k(i, j).s * dx / dy * (j == 1 ? 2.0 : 1.0);
      D.w = k(i, j).w * dy / dx * (i == 1 ? 2.0 : 1.0);

      // Heat flows
      F.n = dx * cp * rho * v(i, j);
      F.e = dy * cp * rho * u(i, j);
      F.s = dx * cp * rho * v(i, j - 1);
      F.w = dy * cp * rho * u(i - 1, j);

      // Compute and store power law coefficients
      fillPowerLaw(T.a(i, j), D, F);
    }
}

void
Problem::uCoefficients()
{
  for (unsigned int i = 1; i < u.Mx; ++i)
    for (unsigned int j = 1; j < u.My; ++j)
    {
      Coefficients D, F;

      // Width of the cell
      const double W = dx * (i == 1 || i == u.Mx - 1 ? 1.5 : 1.0);
      // North/south distances to pressure nodes
      const double dy_pn = dy * (j == u.My - 1 ? 0.5 : 1.0);
      const double dy_ps = dy * (j == 1 ? 0.5 : 1.0);

      // Diffusion coefficients
      D.n = mu_u(i, j).n * W / dy_pn;
      D.e = mu_u(i, j).e * dy / dx;
      D.s = mu_u(i, j).s * W / dy_ps;
      D.w = mu_u(i, j).w * dy / dx;

      // East and west flows
      F.e = rho * dy * (i == u.Mx - 1 ? u(u.Mx, j) : 0.5 * (u(i + 1, j) + u(i, j)));
      F.w = rho * dy * (i == 1 ? u(0, j) : 0.5 * (u(i - 1, j) + u(i, j)));

      // North and south flows
      if (i == 1) // Left boundary
      {
        F.n = rho * W * (v(0, j) + 3.0 * v(1, j) + 2.0 * v(2, j)) / 6.0;
        F.s = rho * W * (v(0, j - 1) + 3.0 * v(1, j - 1) + 2.0 * v(2, j - 1)) / 6.0;
      }
      else if (i == u.Mx - 1) // Right boundary
      {
        F.n = rho * W * (2.0 * v(i, j) + 3.0 * v(i + 1, j) + v(i + 2, j)) / 6.0;
        F.s = rho * W * (2.0 * v(i, j - 1) + 3.0 * v(i + 1, j - 1) + v(i + 2, j - 1)) / 6.0;
      }
      else // Interior (not left or right boundary)
      {
        F.n = rho * W * 0.5 * (v(i, j) + v(i + 1, j));
        F.s = rho * W * 0.5 * (v(i, j - 1) + v(i + 1, j - 1));
      }

      // Pressure RHS
      const double b = dy * (p(i, j) - p(i + 1, j));

      // Compute and store power law coefficients
      fillPowerLaw(u.a(i, j), D, F, b);

      // Explicitly set outflow condition
      if (i == u.Mx - 1)
      {
        u.a(i, j).p -= u.a(i, j).e;
        u.a(i, j).e = 0;
      }
    }
}

void
Problem::vCoefficients()
{
  for (unsigned int i = 1; i < v.Mx; ++i)
    for (unsigned int j = 1; j < v.My; ++j)
    {
      Coefficients D, F;

      // Height of the cell
      const double H = dy * (j == 1 || j == v.My - 1 ? 1.5 : 1.0);
      // East/west distances to pressure nodes
      const double dx_pe = dx * (i == v.Mx - 1 ? 0.5 : 1.0);
      const double dx_pw = dx * (i == 1 ? 0.5 : 1.0);

      // Diffusion coefficient
      D.n = mu_v(i, j).n * dx / dy;
      D.e = mu_v(i, j).e * H / dx_pe;
      D.s = mu_v(i, j).s * dx / dy;
      D.w = mu_v(i, j).w * H / dx_pw;

      // North and east flows
      F.n = rho * dx * (j == v.My - 1 ? v(i, v.My) : 0.5 * (v(i, j + 1) + v(i, j)));
      F.s = rho * dx * (j == 1 ? v(i, 0) : 0.5 * (v(i, j - 1) + v(i, j)));
      // East and west flows
      if (j == 1) // Bottom boundary
      {
        F.e = rho * H * (u(i, 0) + 3.0 * u(i, 1) + 2.0 * u(i, 2)) / 6.0;
        F.w = rho * H * (u(i - 1, 0) + 3.0 * u(i - 1, 1) + 2.0 * u(i - 1, 2)) / 6.0;
      }
      else if (j == v.My - 1) // Top boundary
      {
        F.e = rho * H * (2.0 * u(i, j) + 3.0 * u(i, j + 1) + u(i, j + 2)) / 6.0;
        F.w = rho * H * (2.0 * u(i - 1, j) + 3.0 * u(i - 1, j + 1) + u(i - 1, j + 2)) / 6.0;
      }
      else // Interior (not top or bottom boundary)
      {
        F.e = rho * H * 0.5 * (u(i, j) + u(i, j + 1));
        F.w = rho * H * 0.5 * (u(i - 1, j) + u(i - 1, j + 1));
      }

      // Pressure RHS
      const double b = dx * (p(i, j) - p(i, j + 1));

      // Compute and store power law coefficients
      fillPowerLaw(v.a(i, j), D, F, b);
    }
}

void
Problem::fillPowerLaw(Coefficients & a,
                      const Coefficients & D,
                      const Coefficients & F,
                      const double & b)
{
  a.n = D.n * fmax(0, pow5(1 - 0.1 * fabs(F.n / D.n))) + fmax(-F.n, 0);
  a.e = D.e * fmax(0, pow5(1 - 0.1 * fabs(F.e / D.e))) + fmax(-F.e, 0);
  a.s = D.s * fmax(0, pow5(1 - 0.1 * fabs(F.s / D.s))) + fmax(F.s, 0);
  a.w = D.w * fmax(0, pow5(1 - 0.1 * fabs(F.w / D.w))) + fmax(F.w, 0);
  a.p = a.n + a.e + a.s + a.w;
  a.b = b;
}

} // namespace Flow2D
