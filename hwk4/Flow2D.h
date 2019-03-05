#ifndef Flow2D_H
#define Flow2D_H

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>

#include "Matrix.h"
#include "TriDiagonal.h"

template <typename T>
void
saveCSV(const std::vector<T> & v, std::string filename)
{
  std::ofstream f;
  f.open(filename);
  for (unsigned int i = 0; i < v.size(); ++i)
    f << std::scientific << v[i] << std::endl;
  f.close();
}

template <typename T>
void
print(const std::vector<T> & v, const unsigned int pr = 6)
{
  for (unsigned int i = 0; i < v.size(); ++i)
    std::cout << std::showpos << std::scientific << std::setprecision(pr) << v[i] << " ";
  std::cout << std::endl;
}

struct BoundaryCondition
{
  BoundaryCondition(double top = 0, double right = 0, double bottom = 0, double left = 0)
    : top(top), right(right), bottom(bottom), left(left)
  {
  }
  const double top, right, bottom, left;
};

struct Coefficients
{
  double p, n, e, s, w, b = 0;
  void print(const unsigned int pr = 6)
  {
    std::cout << std::setprecision(pr) << std::scientific << "n = " << n << ", e = " << e
              << ", s = " << s << ", w = " << w << ", p = " << p << ", b = " << b << std::endl;
  }
};

struct MatrixCoefficients
{
  MatrixCoefficients(const unsigned int Nx, const unsigned int Ny) : vals(Nx, Ny), Nx(Nx), Ny(Ny) {}
  Coefficients & operator()(unsigned int i, unsigned int j) { return vals(i, j); }
  Matrix<Coefficients> vals;
  const unsigned int Nx, Ny;
  void print(const std::string prefix = "", const unsigned int pr = 6)
  {
    for (unsigned int i = 1; i < Nx - 1; ++i)
      for (unsigned int j = 1; j < Ny - 1; ++j)
      {
        std::cout << prefix << "(" << i << ", " << j << "): ";
        vals(i, j).print(pr);
      }
  }
};

enum Equations
{
  u,
  v,
  pc
};

static std::string EquationName(Equations eq) {
  switch (eq)
  {
    case Equations::u:
      return "u";
    case Equations::v:
      return "v";
    case Equations::pc:
      return "pc";
  }
}

struct System
{
  System(const Equations eq,
         const unsigned int Nx,
         const unsigned int Ny,
         BoundaryCondition bc = BoundaryCondition())
    : eq(eq),
      name(EquationName(eq)),
      Nx(Nx),
      Ny(Ny),
      Mx(Nx - 1),
      My(Ny - 1),
      a(Nx, Ny),
      phi(Nx, Ny),
      Ax(Nx - 2),
      Ay(Ny - 2),
      bx(Nx - 2),
      by(Ny - 2)
  {
    // Apply initial boundary conditions
    phi.setColumn(0, bc.left);
    phi.setColumn(Mx, bc.right);
    phi.setRow(0, bc.bottom);
    phi.setRow(My, bc.top);
  }

  const double & operator()(unsigned int i, unsigned int j) const { return phi(i, j); }
  double & operator()(unsigned int i, unsigned int j) { return phi(i, j); }
  void reset() { phi = 0; }
  void print(std::string prefix, bool newline = false) { phi.print(prefix, newline); }

  const Equations eq;
  const std::string name;
  const unsigned int Nx, Ny;
  const unsigned int Mx, My;
  MatrixCoefficients a;
  Matrix<double> phi;
  TriDiagonal<double> Ax, Ay;
  std::vector<double> bx, by;
};

/**
 * Solves a 2D heat conduction problem with dirichlet conditions on the top,
 * left, bottom and with a zero-flux condition on the right with Nx x Ny
 * internal control volumes.
 */
class Flow2D
{
public:
  Flow2D(const unsigned int Nx,
         const unsigned int Ny,
         const double Lx,
         const double Ly,
         const BoundaryCondition u_bc,
         const BoundaryCondition v_bc,
         const double rho,
         const double mu,
         const bool loud = false,
         const unsigned int max_its = 1000,
         const double tol = 1.0e-6);

  // Primary runner
  void run();

  // See if this is solved/converged
  bool isConverged() { return converged; }

  // Get the residuals and number of iterations
  const std::vector<double> & getResiduals() const { return residuals; }
  unsigned int getNumIterations() { return residuals.size(); }

private:
  // Flow2D_correct.cpp
  void correct();
  void pCorrect();
  void uCorrect();
  void vCorrect();

  // Flow2D_coefficients.cpp
  void fillCoefficients(System & sys);
  void pcCoefficients();
  void uCoefficients();
  void vCoefficients();
  void velocityCoefficients(Coefficients & a,
                            const Coefficients & D,
                            const Coefficients & F,
                            const double b);

  // Flow2D_residuals.cpp
  void computeResiduals();
  double pResidual();
  double uResidual();
  double vResidual();

  // Flow2D_solvers.cpp
  void solve(System & sys);
  void solveColumn(const unsigned int i, System & sys);
  void solveRow(const unsigned int j, System & sys);
  void solveVelocities();

protected:
  // Number of pressure CVs
  const unsigned int Nx, Ny;

  // Geometry [m]
  const double Lx, Ly, L_ref, dx, dy;
  // Reference velocity
  const double u_ref;
  // Material properties
  const double rho, mu;

  // Loud (debug)
  const bool loud;

  // Maximum iterations
  const unsigned int max_its;
  // Iteration tolerance
  const double tol;
  // Relaxation coefficients
  const double w_uv, alpha_p;

  // Systems
  System u, v, pc;
  // Pressure solution
  Matrix<double> p;

  // Whether or not we converged
  bool converged = false;
  // Residual for each iteration
  std::vector<double> residuals;
};

#endif /* Flow2D_H */
