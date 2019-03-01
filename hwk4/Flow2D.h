#ifndef Flow2D_H
#define Flow2D_H

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "Matrix.h"
#include "TriDiagonal.h"

template <typename T>
void saveCSV(const std::vector<T> &v, std::string filename) {
  std::ofstream f;
  f.open(filename);
  for (unsigned int i = 0; i < v.size(); ++i)
    f << std::scientific << v[i] << std::endl;
  f.close();
}

struct BoundaryCondition {
  BoundaryCondition(double top, double right, double bottom, double left)
      : top(top), right(right), bottom(bottom), left(left) {}
  double top, right, bottom, left;
};

/**
 * Solves a 2D heat conduction problem with dirichlet conditions on the top,
 * left, bottom and with a zero-flux condition on the right with Nx x Ny
 * internal control volumes.
 */
class Flow2D {
public:
  Flow2D(unsigned int Nx, unsigned int Ny, double Lx, double Ly,
         BoundaryCondition u_BC, BoundaryCondition v_BC, double rho, double k,
         double mu, double C_p, unsigned int max_its = 1000);

  void solve();

  // See if this is solved/converged
  bool converged() {
    return (residuals.size() != 0 && residuals.size() != max_its);
  }

  // Get the residuals and number of iterations
  const std::vector<double> &getResiduals() const { return residuals; }
  unsigned int getNumIterations() { return residuals.size(); }

private:
  double computeResidual() const;

  // Precompute operations
  void precomputeProperties();
  void precomputeColumn(unsigned int col);
  void precomputeRow(unsigned int row);

  // Solve and sweep operations
  void solveColumn(unsigned int col);
  void solveRow(unsigned int row);

protected:
  // Number of interior nodal points in the x and y-dimensions
  const unsigned int Nx, Ny;

  // Geometry [m]
  const double Lx, Ly, dx, dy;
  // Boundary conditions
  const BoundaryCondition u_BC, v_BC;
  // Material properties
  const double rho, k, mu, C_p;
  // Properties stored in matrix form
  Matrix<double> a_p, a_n, a_e, a_s, a_w;

  // Maximum iterations
  const unsigned int max_its;

  // Velocity solutions
  Matrix<double> u, v;
  // Pressure solution
  Matrix<double> P;

  // Matrices for the TDMA solves
  TriDiagonal<double> A_x, A_y;
  // RHS/solution vector for the TDMA solves
  std::vector<double> b_x, b_y;
  // Residual for each iteration
  std::vector<double> residuals;
};

#endif /* Flow2D_H */
