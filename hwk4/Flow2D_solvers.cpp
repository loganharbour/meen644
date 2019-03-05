#include "Flow2D.h"

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

  // Solve from left to right
  for (int i = 1; i < Mx; ++i)
    solveColumn(i, eq);
  // Solve from bottom to top
  for (int j = 1; j < My; ++j)
    solveRow(j, eq);

  if (loud)
    switch (eq)
    {
      case Equations::u:
        u.print("u sweep solution = ", true);
        break;
      case Equations::v:
        v.print("v sweep solution = ", true);
        break;
      case Equations::pc:
        pc.print("pc sweep solution = ", true);
        break;
    }
}

void
Flow2D::solveRow(const unsigned int j, const Equations eq)
{
  if (loud)
    std::cout << "Solving row " << j << std::endl;

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

  if (loud)
  {
    std::cout << "A = " << std::endl;
    A.print();
    std::cout << "b = " << std::endl;
    print(b);
  }

  A.solveTDMA(b);

  if (loud)
  {
    std::cout << "sol = " << std::endl;
    print(b);
    std::cout << std::endl;
  }

  for (unsigned int i = 1; i < M; ++i)
    phi(i, j) = b[i - 1];
}

void
Flow2D::solveColumn(const unsigned int i, const Equations eq)
{
  if (loud)
    std::cout << "Solving column " << i << std::endl;

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

  if (loud)
  {
    std::cout << "A = " << std::endl;
    A.print();
    std::cout << "b = " << std::endl;
    print(b);
  }

  A.solveTDMA(b);

  if (loud)
  {
    std::cout << "sol = " << std::endl;
    print(b);
    std::cout << std::endl;
  }

  for (unsigned int j = 1; j < M; ++j)
    phi(i, j) = b[j - 1];
}
