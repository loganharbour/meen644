#ifndef VARIABLE_H
#define VARIABLE_H

#include "Matrix.h"
#include "TriDiagonal.h"
#include "Vector.h"

namespace Flow2D
{

using namespace std;

// Storage for boundary conditions
struct BoundaryCondition
{
  BoundaryCondition() {}
  BoundaryCondition(const double top, const double right, const double bottom, const double left)
    : top(top), right(right), bottom(bottom), left(left)
  {
  }
  bool nonzero() const { return top != 0 || right != 0 || bottom != 0 || left != 0; }
  double top = 0, right = 0, bottom = 0, left = 0;
};

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
  p
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
  }
}

// General storage structure for primary and auxilary variables
struct Variable
{
  // Constructor for a primary variable
  Variable(const Variables name,
           const unsigned int Nx,
           const unsigned int Ny,
           const double alpha,
           const BoundaryCondition bc = BoundaryCondition())
    : name(name),
      string(VariableString(name)),
      Nx(Nx),
      Ny(Ny),
      Mx(Nx - 1),
      My(Ny - 1),
      w(1 / alpha),
      bc(bc),
      a(Nx, Ny),
      phi(Nx, Ny),
      Ax(Nx - 2),
      Ay(Ny - 2),
      bx(Nx - 2),
      by(Ny - 2)
  {
    // Apply initial boundary conditions
    if (bc.left != 0)
      phi.setColumn(0, bc.left);
    if (bc.right != 0)
      phi.setColumn(Mx, bc.right);
    if (bc.bottom != 0)
      phi.setRow(0, bc.bottom);
    if (bc.top != 0)
      phi.setRow(My, bc.top);
  }

  // Constructor for an auxilary variable (no solver storage)
  Variable(const Variables name, const unsigned int Nx, const unsigned int Ny)
    : name(name), string(VariableString(name)), Nx(Nx), Ny(Ny), Mx(Nx - 1), My(Ny - 1), phi(Nx, Ny)
  {
  }

  // Solution matrix operations
  const double & operator()(const unsigned int i, const unsigned int j) const { return phi(i, j); }
  double & operator()(const unsigned int i, const unsigned int j) { return phi(i, j); }
  void print(const string prefix = "", const bool newline = false, const unsigned int pr = 5) const
  {
    phi.print(prefix, newline, pr);
  }
  void save(const string filename) const { phi.save(filename); }
  void reset() { phi = 0; }

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
  // Maximum variable index that is being solved
  const unsigned int Mx, My;
  // Relaxation coefficient used in solving linear systems
  const double w = 0;
  // Boundary conditions
  const BoundaryCondition bc = BoundaryCondition();
  // Matrix coefficients
  Matrix<Coefficients> a;
  // Variable solution
  Matrix<double> phi;
  // Linear system LHS for both sweep directions
  TriDiagonal<double> Ax, Ay;
  // Linear system RHS for both sweep directions
  Vector<double> bx, by;
};

} // namespace Flow2D
#endif /* VARIABLE_H */
