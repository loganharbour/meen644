#include "Problem.h"

namespace Flow2D
{

void
Problem::solveMain()
{
  ++main_iterations;
  if (debug)
    cout << endl << "Main iteration " << main_iterations << endl << endl;

  solve(u);
  solve(v);
  solve(pc);
}

void
Problem::solveAux()
{
  ++aux_iterations;
  if (debug)
    cout << endl << "Aux iteration " << aux_iterations << endl << endl;

  solve(T);
}

void
Problem::solve(Variable & var)
{
  if (debug)
    cout << "Solving variable " << var.string << endl << endl;

  // Fill the coefficients
  fillCoefficients(var);

  // Solve west to east
  sweepColumns(var);
  // Solve south to north
  sweepRows(var);
  // Solve east to west
  sweepColumns(var, false);

  if (debug)
    var.print(var.string + " sweep solution = ", true);
}

void
Problem::sweepRows(Variable & var, const bool south_north)
{
  if (debug)
    cout << "Sweeping " << var.string << (south_north ? " south to north" : " north to south")
         << endl;

  // Sweep south to north
  if (south_north)
    for (int j = 1; j < var.My; ++j)
      sweepRow(j, var);
  // Sweep north to south
  else
    for (int j = var.My - 1; j > 0; --j)
      sweepRow(j, var);
}

void
Problem::sweepColumns(Variable & var, const bool west_east)
{
  if (debug)
    cout << "Sweeping " << var.string << (west_east ? " east to west" : " west to east") << endl;

  // Sweep west to east
  if (west_east)
    for (int i = 1; i < var.Mx; ++i)
      sweepColumn(i, var);
  // Sweep east to west
  else
    for (int i = var.Mx - 1; i > 0; --i)
      sweepColumn(i, var);
}

void
Problem::sweepColumn(const unsigned int i, Variable & var)
{
  if (debug)
    cout << "Solving " << var.string << " column " << i << endl;

  auto & A = var.Ay;
  auto & b = var.by;

  // Fill for each cell
  for (unsigned int j = 1; j < var.My; ++j)
  {
    const Coefficients & a = var.a(i, j);
    b[j - 1] = a.b + a.w * var(i - 1, j) + a.e * var(i + 1, j);
    if (var.w != 1)
      b[j - 1] += a.p * var(i, j) * (var.w - 1);
    if (j == 1)
    {
      A.setTopRow(a.p * var.w, -a.n);
      if (var.name != Variables::pc)
        b[j - 1] += a.s * var(i, j - 1);
    }
    else if (j == var.My - 1)
    {
      A.setBottomRow(-a.s, a.p * var.w);
      if (var.name != Variables::pc)
        b[j - 1] += a.n * var(i, j + 1);
    }
    else
      A.setMiddleRow(j - 1, -a.s, a.p * var.w, -a.n);
  }

  if (debug)
  {
    A.print("A =");
    b.print("b =");
  }

  // Solve
  A.TDMA(b);

  if (debug)
    b.print("sol =", true);

  // Store solution
  for (unsigned int j = 1; j < var.My; ++j)
    var(i, j) = b[j - 1];
}

void
Problem::sweepRow(const unsigned int j, Variable & var)
{
  if (debug)
    cout << "Solving " << var.string << " row " << j << endl;

  auto & A = var.Ax;
  auto & b = var.bx;

  // Fill for each cell
  for (unsigned int i = 1; i < var.Mx; ++i)
  {
    const Coefficients & a = var.a(i, j);
    b[i - 1] = a.b + a.s * var(i, j - 1) + a.n * var(i, j + 1);
    if (var.w != 1)
      b[i - 1] += a.p * var(i, j) * (var.w - 1);
    if (i == 1)
    {
      A.setTopRow(a.p * var.w, -a.e);
      if (var.name != Variables::pc)
        b[i - 1] += a.w * var(i - 1, j);
    }
    else if (i == var.Mx - 1)
    {
      A.setBottomRow(-a.w, a.p * var.w);
      if (var.name != Variables::pc)
        b[i - 1] += a.e * var(i + 1, j);
    }
    else
      A.setMiddleRow(i - 1, -a.w, a.p * var.w, -a.e);
  }

  if (debug)
  {
    A.print("A =");
    b.print("b =");
  }

  // Solve
  A.TDMA(b);

  if (debug)
    b.print("sol =", true);

  // Store solution
  for (unsigned int i = 1; i < var.Mx; ++i)
    var(i, j) = b[i - 1];
}

} // namespace Flow2D
