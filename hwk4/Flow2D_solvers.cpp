#include "Flow2D.h"

void
Flow2D::solve(System & sys)
{
  unsigned int Mx, My;
  fillCoefficients(sys);

  // Solve from left to right
  for (int i = 1; i < sys.Mx; ++i)
    solveColumn(i, sys);
  // Solve from bottom to top
  for (int j = 1; j < sys.My; ++j)
    solveRow(j, sys);

  if (loud)
    sys.print(sys.name + " sweep solution = ", true);
}

void
Flow2D::solveRow(const unsigned int j, System & sys)
{
  if (loud)
    std::cout << "Solving " << sys.name << " row " << j << std::endl;

  auto & A = sys.Ax;
  auto & b = sys.bx;
  double w = (sys.eq == Equations::pc ? 1 : w_uv);

  // Fill for each cell
  for (unsigned int i = 1; i < sys.Mx; ++i)
  {
    Coefficients & a = sys.a(i, j);
    b[i - 1] = a.b + a.s * sys(i, j - 1) + a.n * sys(i, j + 1) + a.p * sys(i, j) * (w - 1);
    if (i == 1)
    {
      A.setTopRow(a.p * w, -a.e);
      if (sys.eq != Equations::pc)
        b[i - 1] += a.w * sys(i - 1, j);
    }
    else if (i == sys.Mx - 1)
    {
      A.setBottomRow(-a.w, a.p * w);
      if (sys.eq != Equations::pc)
        b[i - 1] += a.e * sys(i + 1, j);
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

  for (unsigned int i = 1; i < sys.Mx; ++i)
    sys(i, j) = b[i - 1];
}

void
Flow2D::solveColumn(const unsigned int i, System & sys)
{
  if (loud)
    std::cout << "Solving " << sys.name << " column " << i << std::endl;

  // Grab the appropriate matrix/vector/coefficients/previous solution
  auto & A = sys.Ay;
  auto & b = sys.by;
  double w = (sys.eq == Equations::pc ? 1 : w_uv);

  // Fill for each cell
  for (unsigned int j = 1; j < sys.My; ++j)
  {
    Coefficients & a = sys.a(i, j);
    b[j - 1] = a.b + a.w * sys(i - 1, j) + a.e * sys(i + 1, j) + a.p * sys(i, j) * (w - 1);
    if (j == 1)
    {
      A.setTopRow(a.p * w, -a.n);
      if (sys.eq != Equations::pc)
        b[j - 1] += a.s * sys(i, j - 1);
    }
    else if (j == sys.My - 1)
    {
      A.setBottomRow(-a.s, a.p * w);
      if (sys.eq != Equations::pc)
        b[j - 1] += a.n * sys(i, j + 1);
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

  for (unsigned int j = 1; j < sys.My; ++j)
    sys(i, j) = b[j - 1];
}
