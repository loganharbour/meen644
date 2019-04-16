#ifndef VARIABLE_H
#define VARIABLE_H

#include "Matrix.h"
#include "TriDiagonal.h"
#include "Vector.h"

#include <functional>

namespace Flow2D
{

using namespace std;

// Storage for coefficients for a single CV
struct Coefficients
{
  double p = 0, n = 0, e = 0, s = 0, w = 0, b = 0;
  void print(const unsigned int pr = 5) const
  {
    cout << setprecision(pr) << scientific << "n = " << n << ", e = " << e << ", s = " << s
         << ", w = " << w << ", p = " << p << ", b = " << b << endl;
  }
};

// Enum for variable types
enum Variables
{
  u,
  v,
  pc,
  p,
  T
};

// Conversion from variable type to its string
static string
VariableString(Variables var)
{
  switch (var)
  {
    case Variables::u:
      return "u";
    case Variables::v:
      return "v";
    case Variables::pc:
      return "pc";
    case Variables::p:
      return "p";
    case Variables::T:
      return "T";
  }
}

// General storage structure for primary and auxilary variables
struct Variable
{
  // Constructor for a primary variable
  Variable(const Variables name,
           const unsigned int Nx,
           const unsigned int Ny,
           const double dx,
           const double dy,
           const double alpha,
           function<double(const vector<double>)> ic = [](const vector<double> p) { return 0; })
    : name(name),
      string(VariableString(name)),
      Nx(Nx),
      Ny(Ny),
      dx(dx),
      dy(dy),
      Mx(Nx - 1),
      My(Ny - 1),
      w(1 / alpha),
      a(Nx, Ny),
      phi(Nx, Ny),
      Ax(Nx - 2),
      Ay(Ny - 2),
      bx(Nx - 2),
      by(Ny - 2)
  {
    for (unsigned int i = 0; i <= Mx; ++i)
    {
      phi(i, 0) = ic(point(i, 0));
      phi(i, My) = ic(point(i, My));
    }
    for (unsigned int j = 0; j <= My; ++j)
    {
      phi(0, j) = ic(point(0, j));
      phi(Mx, j) = ic(point(Mx, j));
    }
  }

  // Constructor for an auxilary variable (no solver storage)
  Variable(const Variables name,
           const unsigned int Nx,
           const unsigned int Ny,
           const double dx,
           const double dy,
           function<double(const vector<double>)> ic = [](const vector<double> p) { return 0; })
    : name(name),
      string(VariableString(name)),
      Nx(Nx),
      Ny(Ny),
      dx(dx),
      dy(dy),
      Mx(Nx - 1),
      My(Ny - 1),
      phi(Nx, Ny)
  {
    for (unsigned int i = 0; i <= Mx; ++i)
    {
      phi(i, 0) = ic(point(i, 0));
      phi(i, My) = ic(point(i, My));
    }
    for (unsigned int j = 0; j <= My; ++j)
    {
      phi(0, j) = ic(point(0, j));
      phi(Mx, j) = ic(point(Mx, j));
    }
  }

  // Solution matrix operations
  void operator=(const double v) { phi = v; }
  const double & operator()(const unsigned int i, const unsigned int j) const { return phi(i, j); }
  double & operator()(const unsigned int i, const unsigned int j) { return phi(i, j); }
  void print(const string prefix = "", const bool newline = false, const unsigned int pr = 5) const
  {
    phi.print(prefix, newline, pr);
  }
  void save(const string filename) const { phi.save(filename); }
  void reset() { phi = 0; }

  // Get the point in space associated with an i, j for this CV
  const vector<double> point(const unsigned int i, const unsigned int j) const
  {
    const double id = (double)i, jd = (double)j;
    vector<double> pt = {dx * (i != 0) * (id - 0.5 * (1 + (i == Mx))),
                         dy * (j != 0) * (jd - 0.5 * (1 + (j == My)))};
    if (name == Variables::u)
      pt[0] = id * dx;
    if (name == Variables::v)
      pt[1] = jd * dy;
    return pt;
  }

  // Coefficient debug
  void printCoefficients(const string prefix = "",
                         const bool newline = false,
                         const unsigned int pr = 5) const
  {
    for (unsigned int i = 1; i < Nx - 1; ++i)
      for (unsigned int j = 1; j < Ny - 1; ++j)
      {
        cout << prefix << "(" << i << ", " << j << "): ";
        a(i, j).print(pr);
      }
    if (newline)
      cout << endl;
  }

  // Variable enum name
  const Variables name;
  // Variable string
  const string string;
  // Variable size
  const unsigned int Nx, Ny;
  // Mesh size
  const double dx, dy;
  // Maximum variable index that is being solved
  const unsigned int Mx, My;
  // Relaxation coefficient used in solving linear systems
  const double w = 0;
  // Matrix coefficients
  Matrix<Coefficients> a;
  // Variable solution
  Matrix<double> phi;
  // Linear system LHS for both sweep directions
  TriDiagonal Ax, Ay;
  // Linear system RHS for both sweep directions
  Vector<double> bx, by;
}; // namespace Flow2D

} // namespace Flow2D
#endif /* VARIABLE_H */
