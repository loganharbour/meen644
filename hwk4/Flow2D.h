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

struct BoundaryCondition
{
  BoundaryCondition(double top, double right, double bottom, double left)
    : top(top), right(right), bottom(bottom), left(left)
  {
  }
  double top, right, bottom, left;
};

struct Coefficients
{
  double p, n, e, s, w, b = 0;
};

struct MatrixCoefficients
{
  MatrixCoefficients(unsigned int Nx, unsigned int Ny) : vals(Nx, Ny) {}
  Coefficients & operator()(unsigned int i, unsigned int j) { return vals(i, j); }
  Matrix<Coefficients> vals;
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
  Flow2D(unsigned int Nx,
         unsigned int Ny,
         double Lx,
         double Ly,
         BoundaryCondition u_BC,
         BoundaryCondition v_BC,
         double rho,
         double mu,
         unsigned int max_its = 1000);

  void solve();
  // See if this is solved/converged
  bool converged() { return (residuals.size() != 0 && residuals.size() != max_its); }

  // Get the residuals and number of iterations
  const std::vector<double> & getResiduals() const { return residuals; }
  unsigned int getNumIterations() { return residuals.size(); }

private:
  void correct();
  void fillInitialValues();

  void uCoefficients();
  void vCoefficients();
  void pcCoefficients();
  void velocityCoefficients(Coefficients & a,
                            const Coefficients & D,
                            const Coefficients & F,
                            const double b);

  double uResidual();
  double vResidual();
  double pResidual();

  void solve(const Equations eq);

  // Solve and sweep operations
  void solveVelocities();
  void solveColumn(const unsigned int i, const Equations equation);
  void solveRow(const unsigned int j, const Equations equation);

protected:
  // Number of pressure CVs
  const unsigned int Nx, Ny;
  // Maximum nodal values
  const unsigned int Mx_u, My_u, Mx_v, My_v, Mx_p, My_p;

  // Geometry [m]
  const double Lx, Ly, dx, dy;
  // Boundary conditions
  const BoundaryCondition u_BC, v_BC;
  // Material properties
  const double rho, mu;
  // Coefficient matrices
  MatrixCoefficients a_u, a_v, a_pc;

  // Maximum iterations
  const unsigned int max_its;
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

  // Residual for each iteration
  std::vector<double> residuals;
};

#endif /* Flow2D_H */
