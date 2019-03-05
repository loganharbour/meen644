#include "Problem.h"

namespace Flow2D
{

void
Problem::computeResiduals()
{
  double Ru = velocityResidual(u);
  double Rv = velocityResidual(v);
  double Rp = pResidual();

  // Check for convergence
  if (Ru < tol && Rv < tol && Rp < tol)
    converged = true;

  // Print residuals
  if (converged)
    cout << "Converged in " << iterations << " iterations: ";
  if (converged || loud) {
    cout << setprecision(2) << scientific;
    cout << "u = " << Ru;
    cout << ", v = " << Rv;
    cout << ", p = " << Rp << endl;
  }
}

double
Problem::pResidual()
{
  double numer = 0;
  for (unsigned int i = 1; i < pc.Mx; ++i)
    for (unsigned int j = 1; j < pc.My; ++j)
      numer += abs(dy * (u(i - 1, j) - u(i, j)) + dx * (v(i, j - 1) - v(i, j)));
  return numer / (u_ref * L_ref);
}

double
Problem::velocityResidual(const Variable & var)
{
  double numer, numer_temp, denom = 0;
  for (unsigned int i = 1; i < var.Mx; ++i)
    for (unsigned int j = 1; j < var.My; ++j)
    {
      const Coefficients & a = var.a(i, j);
      numer_temp = a.p * var(i, j);
      denom += abs(numer_temp);
      numer_temp -= a.n * var(i, j + 1);
      numer_temp -= a.e * var(i + 1, j);
      numer_temp -= a.s * var(i, j - 1);
      numer_temp -= a.w * var(i - 1, j);
      numer_temp -= a.b;
      numer += abs(numer_temp);
    }
  return numer / denom;
}

}
