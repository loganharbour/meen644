#ifndef CONDUCTION2D_H
#define CONDUCTION2D_H

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

#include "Matrix.h"
#include "TriDiagonal.h"

class Conduction2D {
public:
  Conduction2D(unsigned int Nx, unsigned int Ny, double alpha, double Lx = 0.5,
               double Ly = 0.5, double k = 386.0, double T_L = 50.0,
               double T_T = 100.0, double T_B = 50.0, double tol = 1e-5,
               unsigned int max_its = 1000);

  void run();
  void save(std::string filename);

  unsigned int getNx() { return Nx; }
  unsigned int getNy() { return Ny; }
  Matrix & getT() { return T; }

private:
  void precomputeProperties();
  void precompute();
  void precomputeColumn(unsigned int col);
  void precomputeRow(unsigned int row);
  void solveColumn(unsigned int col);
  void solveRow(unsigned int row);
  void sweep();
  double computeResidual();

protected:
  // Number of interior nodal points in the x and y-dimensions
  const unsigned int Nx, Ny;

  // Geometry
  const double Lx, Ly, dx, dy;
  // Material properties
  const double k;
  // Boundary conditions
  const double T_L, T_T, T_B;
  // Properties stored in matrix form
  Matrix a_p, a_n, a_e, a_s, a_w;

  // Relaxation coefficient
  const double w_inv;
  // Iteration tolerance
  const double tol;
  // Maximum iterations
  const unsigned int max_its;

  // Temperature solution
  Matrix T;

  // Precomputed matrices for the TDMA solves
  std::vector<TriDiagonal> pre_A_x, pre_A_y;
  // Precomputed RHS for the TDMA solves
  std::vector<std::vector<double> > pre_b_x, pre_b_y;
  // Matrices for the TDMA solves
  TriDiagonal A_x, A_y;
  // RHS/solution vector for the TDMA solves
  std::vector<double> b_x, b_y;
};

#endif /* CONDUCTION2D_H */
