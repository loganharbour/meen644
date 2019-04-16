#include "Problem.h"

namespace Flow2D
{

void
Problem::computeMainResiduals()
{
  const double Rp = pResidual();
  const double Ru = velocityResidual(u);
  const double Rv = velocityResidual(v);

  if (Ru < tol && Rv < tol && Rp < tol)
  {
    if (debug)
      cout << "Main variables converged in " << main_iterations << " iterations" << endl;
    main_converged = true;
  }
}

void
Problem::computeAuxResiduals()
{
  double RT = TResidual();

  // Still not converged
  if (RT > tol)
    return;

  // Converged, finish up
  converged = true;

  // Print the result
  const double Rp = pResidual();
  const double Ru = velocityResidual(u);
  const double Rv = velocityResidual(v);
  cout << "Converged in " << noshowpos << fixed << setprecision(3)
       << 1.0 * (clock() - start) / CLOCKS_PER_SEC << " sec in " << main_iterations << " main and "
       << aux_iterations << " aux iterations: ";
  cout << noshowpos << setprecision(1) << scientific;
  cout << "p = " << Rp;
  cout << ", u = " << Ru;
  cout << ", v = " << Rv;
  cout << ", T = " << RT << endl;
}

double
Problem::pResidual() const
{
  double numer = 0;
  for (unsigned int i = 1; i < pc.Mx; ++i)
    for (unsigned int j = 1; j < pc.My; ++j)
      numer += abs(dy * (u(i - 1, j) - u(i, j)) + dx * (v(i, j - 1) - v(i, j)));
  return numer / (u_ref * L_ref);
}

double
Problem::TResidual() const
{
  double numer, numer_temp, denom = 0;
  for (unsigned int i = 1; i < T.Mx; ++i)
    for (unsigned int j = 1; j < T.My; ++j)
    {
      const Coefficients & a = T.a(i, j);
      numer_temp = a.p * T(i, j);
      denom += abs(numer_temp);
      numer_temp -= a.n * T(i, j + 1) + a.e * T(i + 1, j);
      numer_temp -= a.s * T(i, j - 1) + a.w * T(i - 1, j);
      numer += abs(numer_temp);
    }
  return numer / denom;
}

double
Problem::velocityResidual(const Variable & var) const
{
  double numer, numer_temp, denom = 0;
  for (unsigned int i = 1; i < var.Mx; ++i)
    for (unsigned int j = 1; j < var.My; ++j)
    {
      const Coefficients & a = var.a(i, j);
      numer_temp = a.p * var(i, j);
      denom += abs(numer_temp);
      numer_temp -= a.n * var(i, j + 1) + a.e * var(i + 1, j);
      numer_temp -= a.s * var(i, j - 1) + a.w * var(i - 1, j) + a.b;
      numer += abs(numer_temp);
    }
  return numer / denom;
}

} // namespace Flow2D
