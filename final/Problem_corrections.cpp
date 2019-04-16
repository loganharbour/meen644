#include "Problem.h"

namespace Flow2D
{

void
Problem::correctMain()
{
  uCorrect();
  vCorrect();
  pCorrect();
  pBCCorrect();
  uBCCorrect();
}

void
Problem::correctAux()
{
  TBCCorrect();
}

void
Problem::pCorrect()
{
  for (unsigned int i = 1; i < pc.Mx; ++i)
    for (unsigned int j = 1; j < pc.My; ++j)
      p(i, j) += alpha_p * pc(i, j);

  // Set pressure correction back to zero
  pc.reset();

  if (debug)
    p.print("p corrected = ", true);
}

void
Problem::pBCCorrect()
{
  // Apply the edge values as velocity is set
  for (unsigned int i = 0; i <= pc.Mx; ++i)
  {
    p(i, 0) = p(i, 1);
    p(i, pc.My) = p(i, pc.My - 1);
  }
  for (unsigned int j = 0; j <= pc.My; ++j)
  {
    p(0, j) = p(1, j);
    p(pc.Mx, j) = p(pc.Mx - 1, j);
  }

  if (debug)
    p.print("p boundary condition corrected = ", true);
}
void
Problem::TBCCorrect()
{
  for (unsigned int j = 0; j <= T.My; ++j)
    T(T.Mx, j) = (3 * T(T.Mx - 1, j) - T(T.Mx - 2, j)) / 2;

  for (unsigned int i = 1; i < T.Mx; ++i)
  {
    T(i, 0) = T(i, 1) + q(T.point(i, 0)) * dy / (2 * k(i, 1).p);
    T(i, T.My) = T(i, T.My - 1) + q(T.point(i, T.My)) * dy / (2 * k(i, T.My - 1).p);
  }
}

void
Problem::uCorrect()
{
  for (unsigned int i = 1; i < u.Mx; ++i)
    for (unsigned int j = 1; j < u.My; ++j)
      u(i, j) += dy * (pc(i, j) - pc(i + 1, j)) / u.a(i, j).p;

  if (debug)
    u.print("u corrected = ", true);
}

void
Problem::uBCCorrect()
{
  double m_in = 0, m_out = 0;
  for (unsigned int j = 1; j < u.My; ++j)
  {
    m_in += rho * dy * u(0, j);
    m_out += rho * dy * u(u.Mx - 1, j);
  }
  for (unsigned int j = 0; j <= u.My; ++j)
    u(u.Mx, j) = m_in * u(u.Mx - 1, j) / m_out;

  if (debug)
    u.print("u boundary condition corrected = ", true);
}

void
Problem::vCorrect()
{
  for (unsigned int i = 1; i < v.Mx; ++i)
    for (unsigned int j = 1; j < v.My; ++j)
      v(i, j) += dx * (pc(i, j) - pc(i, j + 1)) / v.a(i, j).p;

  if (debug)
    v.print("v corrected = ", true);
}

} // namespace Flow2D
