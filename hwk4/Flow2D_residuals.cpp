#include "Flow2D.h"

void
Flow2D::computeResiduals()
{
  double Ru = uResidual();
  double Rv = vResidual();
  double Rp = pResidual();

  std::cout << "Ru = " << std::setprecision(4) << std::scientific << Ru
            << ", Rv = " << std::setprecision(4) << std::scientific << Rv
            << ", Rp = " << std::setprecision(4) << std::scientific << Rv << std::endl;

  if (Ru < tol && Rv < tol && Rp < tol)
    converged = true;
}

double
Flow2D::pResidual()
{
  double numer, numer_temp = 0;
  for (unsigned int i = 1; i < Mx_p; ++i)
    for (unsigned int j = 1; j < My_p; ++j)
    {
      numer_temp = dy * (u(i - 1, j) - u(i, j)) + dx * (v(i, j - 1) - v(i, j));
      numer += std::abs(rho * numer_temp);
    }
  return numer / (rho * u_ref * L_ref);
}

double
Flow2D::uResidual()
{
  double numer, numer_temp, denom = 0;
  for (unsigned int i = 1; i < Mx_u; ++i)
    for (unsigned int j = 1; j < My_u; ++j)
    {
      Coefficients & a = a_u(i, j);
      numer_temp = a.p * u(i, j);
      denom += std::abs(numer_temp);
      numer_temp -= a.n * u(i, j + 1);
      numer_temp -= a.e * u(i + 1, j);
      numer_temp -= a.s * u(i, j - 1);
      numer_temp -= a.w * u(i - 1, j);
      numer_temp -= a.b;
      numer_temp -= dy * (p(i, j) - p(i + 1, j));
      numer += std::abs(numer_temp);
    }
  return numer / denom;
}

double
Flow2D::vResidual()
{
  double numer, numer_temp, denom = 0;
  for (unsigned int i = 1; i < Mx_v; ++i)
    for (unsigned int j = 1; j < My_v; ++j)
    {
      Coefficients & a = a_v(i, j);
      numer_temp = a.p * v(i, j);
      denom += std::abs(numer_temp);
      numer_temp -= a.n * v(i, j + 1);
      numer_temp -= a.e * v(i + 1, j);
      numer_temp -= a.s * v(i, j - 1);
      numer_temp -= a.w * v(i - 1, j);
      numer_temp -= a.b;
      numer_temp -= dx * (p(i, j) - p(i, j + 1));
      numer += std::abs(numer_temp);
    }
  return numer / denom;
}
