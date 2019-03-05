#ifndef Flow2D_H
#define Flow2D_H

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

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
  BoundaryCondition(double top, double right, double bottom, double left)
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
         const BoundaryCondition u_BC,
         const BoundaryCondition v_BC,
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
  void fillInitialValues();

  // Flow2D_correct.cpp
  void correct();
  void pCorrect();
  void uCorrect();
  void vCorrect();

  // Flow2D_coefficients.cpp
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
  void solve(const Equations eq);
  void solveColumn(const unsigned int i, const Equations equation);
  void solveRow(const unsigned int j, const Equations equation);
  void solveVelocities();

protected:
  // Number of pressure CVs
  const unsigned int Nx, Ny;
  // Maximum nodal values
  const unsigned int Mx_u, My_u, Mx_v, My_v, Mx_p, My_p;

  // Geometry [m]
  const double Lx, Ly, L_ref, dx, dy;
  // Boundary conditions
  const BoundaryCondition u_BC, v_BC;
  // Reference velocity
  const double u_ref;
  // Material properties
  const double rho, mu;
  // Coefficient matrices
  MatrixCoefficients a_u, a_v, a_pc;

  // Loud (debug)
  const bool loud;

  // Maximum iterations
  const unsigned int max_its;
  // Iteration tolerance
  const double tol;
  // Relaxation coefficients
  const double w_uv, alpha_p;

  // Velocity solutions
  Matrix<double> u, v;
  // Pressure solution
  Matrix<double> p;
  // Pressure corrector solution
  Matrix<double> pc;

  // Matrices and vectors for sweeping
  TriDiagonal<double> Ax_u, Ay_u, Ax_v, Ay_v, Ax_pc, Ay_pc;
  std::vector<double> bx_u, by_u, bx_v, by_v, bx_pc, by_pc;

  // Whether or not we converged
  bool converged = false;
  // Residual for each iteration
  std::vector<double> residuals;
};

#endif /* Flow2D_H */
