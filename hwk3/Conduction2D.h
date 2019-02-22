#ifndef CONDUCTION2D_H
#define CONDUCTION2D_H

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "Matrix.h"
#include "TriDiagonal.h"

template <typename T> void saveCSV(const std::vector<T> &v, std::string filename) {
  std::ofstream f;
  f.open(filename);
  for (unsigned int i = 0; i < v.size(); ++i)
    f << std::scientific << v[i] << std::endl;
  f.close();
}

/**
 * Solves a 2D heat conduction problem with dirichlet conditions on the top,
 * left, bottom and with a zero-flux condition on the right with Nx x Ny
 * internal control volumes.
 */
class Conduction2D {
public:
  Conduction2D(unsigned int Nx, unsigned int Ny, double alpha, double Lx = 0.5,
               double Ly = 0.5, double k = 386.0, double T_L = 50.0,
               double T_T = 100.0, double T_B = 50.0, double tol = 1e-5,
               unsigned int max_its = 2000);

  void solve();

  // See if this is solved/converged
  bool converged() {
    return (residuals.size() != 0 && residuals.size() != max_its);
  }

  // Get the solution at the (i, j) internal node
  const double getT(unsigned int i, unsigned int j) const { return T(i, j); }

  // Get the residuals and number of iterations
  const std::vector<double> &getResiduals() const { return residuals; }
  unsigned int getNumIterations() { return residuals.size(); }

  // Save the solution
  void saveT(std::string filename) const { T.save(filename); }

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
  // Heat conduction coefficient [W / m k]
  const double k;
  // Dirichlet oundary conditions (left, top, bottom) [C]
  const double T_L, T_T, T_B;
  // Properties stored in matrix form
  Matrix<double> a_p, a_n, a_e, a_s, a_w;

  // Inverse of the relaxation coefficient
  const double w_inv;
  // Iteration tolerance
  const double tol;
  // Maximum iterations
  const unsigned int max_its;

  // Temperature solution
  Matrix<double> T;

  // Precomputed matrices for the TDMA solves
  std::vector<TriDiagonal<double>> pre_A_x, pre_A_y;
  // Precomputed RHS for the TDMA solves
  std::vector<std::vector<double>> pre_b_x, pre_b_y;
  // Matrices for the TDMA solves
  TriDiagonal<double> A_x, A_y;
  // RHS/solution vector for the TDMA solves
  std::vector<double> b_x, b_y;
  // Residual for each iteration
  std::vector<double> residuals;
};

#endif /* CONDUCTION2D_H */
