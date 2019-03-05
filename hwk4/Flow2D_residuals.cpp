#include "Flow2D.h"

void
Flow2D::computeResiduals()
{
  double Ru = uResidual();
  double Rv = vResidual();
  double Rp = pResidual();

  std::cout << std::setprecision(2) << std::scientific;
  std::cout << "u = " << Ru;
  std::cout << ", v = " << Rv;
  std::cout << ", p = " << Rp << std::endl;

  if (Ru < tol && Rv < tol && Rp < tol)
    converged = true;
}

double
Flow2D::pResidual()
{
  double numer = 0;
  for (unsigned int i = 1; i < pc.Mx; ++i)
    for (unsigned int j = 1; j < pc.My; ++j)
      numer += std::abs(dy * (u(i - 1, j) - u(i, j)) + dx * (v(i, j - 1) - v(i, j)));
  return numer / (u_ref * L_ref);
}

double
Flow2D::uResidual()
{
  double numer, numer_temp, denom = 0;
  for (unsigned int i = 1; i < u.Mx; ++i)
    for (unsigned int j = 1; j < u.My; ++j)
    {
      Coefficients & a = u.a(i, j);
      numer_temp = a.p * u(i, j);
      denom += std::abs(numer_temp);
      numer_temp -= a.n * u(i, j + 1);
      numer_temp -= a.e * u(i + 1, j);
      numer_temp -= a.s * u(i, j - 1);
      numer_temp -= a.w * u(i - 1, j);
      numer_temp -= a.b;
      numer += std::abs(numer_temp);
    }
  return numer / denom;
}

double
Flow2D::vResidual()
{
  double numer, numer_temp, denom = 0;
  for (unsigned int i = 1; i < v.Mx; ++i)
    for (unsigned int j = 1; j < v.My; ++j)
    {
      Coefficients & a = v.a(i, j);
      numer_temp = a.p * v(i, j);
      denom += std::abs(numer_temp);
      numer_temp -= a.n * v(i, j + 1);
      numer_temp -= a.e * v(i + 1, j);
      numer_temp -= a.s * v(i, j - 1);
      numer_temp -= a.w * v(i - 1, j);
      numer_temp -= a.b;
      numer += std::abs(numer_temp);
    }
  return numer / denom;
}
