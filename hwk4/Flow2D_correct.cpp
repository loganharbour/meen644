#include "Flow2D.h"

void
Flow2D::correct()
{
  uCorrect();
  vCorrect();
  pCorrect();
}

void
Flow2D::pCorrect()
{
  for (unsigned int i = 1; i < Mx_p; ++i)
    for (unsigned int j = 1; j < My_p; ++j)
      p(i, j) += alpha_p * pc(i, j);
  pc = 0;

  if (loud)
    p.print("p corrected = ", true);
}

void
Flow2D::uCorrect()
{
  for (unsigned int i = 1; i < Mx_u; ++i)
    for (unsigned int j = 1; j < My_u; ++j)
      u(i, j) += dy * (pc(i, j) - pc(i + 1, j)) / a_u(i, j).p;
  if (loud)
    u.print("u corrected = ", true);
}

void
Flow2D::vCorrect()
{
  for (unsigned int i = 1; i < Mx_v; ++i)
    for (unsigned int j = 1; j < My_v; ++j)
      v(i, j) += dx * (pc(i, j) - pc(i, j + 1)) / a_v(i, j).p;

  if (loud)
    v.print("v corrected = ", true);
}

void
Flow2D::pBoundaryCorrect()
{
  for (unsigned int i = 0; i <= Mx_p; ++i)
  {
    p(i, 0) = p(i, 1);
    p(i, My_p) = p(i, My_p - 1);
  }
  for (unsigned int j = 0; j <= My_p; ++j)
  {
    p(0, j) = p(1, j);
    p(Mx_p, j) = p(Mx_p - 1, j);
  }
}
