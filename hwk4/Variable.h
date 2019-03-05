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
  BoundaryCondition(double top = 0, double right = 0, double bottom = 0, double left = 0)
    : top(top), right(right), bottom(bottom), left(left)
  {
  }
  bool nonzero() const { return top != 0 || right != 0 || bottom != 0 || left != 0; }
  double top, right, bottom, left;
};

// Storage for coefficients for a single CV
struct Coefficients
{
  double p, n, e, s, w, b = 0;
  void print(const unsigned int pr = 6)
  {
    cout << setprecision(pr) << scientific << "n = " << n << ", e = " << e << ", s = " << s
         << ", w = " << w << ", p = " << p << ", b = " << b << endl;
  }
};

// Storage for coefficients for every CV
struct MatrixCoefficients
{
  MatrixCoefficients() {}
  MatrixCoefficients(const unsigned int Nx, const unsigned int Ny) : vals(Nx, Ny), Nx(Nx), Ny(Ny) {}
  Coefficients & operator()(unsigned int i, unsigned int j) { return vals(i, j); }
  const Coefficients & operator()(unsigned int i, unsigned int j) const { return vals(i, j); }
  Matrix<Coefficients> vals;
  const unsigned int Nx = 0, Ny = 0;
  void print(const string prefix = "", const unsigned int pr = 6)
  {
    for (unsigned int i = 1; i < Nx - 1; ++i)
      for (unsigned int j = 1; j < Ny - 1; ++j)
      {
        cout << prefix << "(" << i << ", " << j << "): ";
        vals(i, j).print(pr);
      }
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

struct Variable
{
  // Constructor for a solved variable
  Variable(const Variables name,
           const unsigned int Nx,
           const unsigned int Ny,
           const double alpha,
           BoundaryCondition bc = BoundaryCondition())
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

  // Constructor for an aux variable
  Variable(const Variables name,
           const unsigned int Nx,
           const unsigned int Ny)
    : name(name), string(VariableString(name)), Nx(Nx), Ny(Ny), Mx(Nx - 1), My(Ny - 1), phi(Nx, Ny)
  {
  }

  // Matrix operations
  const double & operator()(unsigned int i, unsigned int j) const { return phi(i, j); }
  double & operator()(unsigned int i, unsigned int j) { return phi(i, j); }
  void print(const string prefix = "", const bool newline = false) const
  {
    phi.print(prefix, newline);
  }
  void save(const string filename) const { phi.save(filename); }
  void reset() { phi = 0; }

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
  MatrixCoefficients a;
  // Variable solution
  Matrix<double> phi;
  // Linear system LHS
  TriDiagonal<double> Ax, Ay;
  // Linear system RHS
  Vector<double> bx, by;
};

}
#endif /* VARIABLE_H */
