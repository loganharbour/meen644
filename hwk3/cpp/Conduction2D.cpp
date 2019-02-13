#include "Conduction2D.h"

Conduction2D::Conduction2D(unsigned int Nx, unsigned int Ny, double w,
                           double Lx, double Ly, double k, double T_L,
                           double T_T, double T_B, double tol,
                           unsigned int max_its)
    : // Interior nodal points
      Nx(Nx), Ny(Ny),
      // Domain size
      Lx(Lx), Ly(Ly), dx(Lx / Nx), dy(Ly / Ny), // [m]
      // Material properties
      k(k), // [W / m - K]
      // Boundary conditions
      T_L(T_L), T_T(T_T), T_B(T_B),
      // Material properties in matrix form
      a_p(Nx, Ny), a_n(Nx, Ny), a_e(Nx, Ny), a_s(Nx, Ny), a_w(Nx, Ny),
      // Relaxation coefficient and iteration tolerance
      w_inv(1.0 / w), tol(tol), max_its(max_its),
      // Initialize matrices and vectors
      // T(Nx, Ny, (T_L + T_T + T_B) / 3.0)
      T(Nx, Ny, 0), pre_A_x(Ny, Nx), pre_A_y(Ny, Nx),
      pre_b_x(Ny), pre_b_y(Nx), A_x(Nx), A_y(Ny), b_x(Nx), b_y(Ny) {}

void Conduction2D::run() {
  precompute();
  for (unsigned int l = 1; l <= max_its; ++l) {
    sweep();

    double R = computeResidual();
    if (R < tol)
    {
      std::cout << "Converged with " << l << " iterations" << std::endl;
      break;
    }
    if (l == max_its)
      std::cout << "Failed to converge" << std::endl;
  }
}

void Conduction2D::save(std::string filename) { T.save(filename); }

void Conduction2D::precomputeProperties() {
  // Set all neighbors to the default at first
  a_n = k * dx / dy;
  a_e = k * dy / dx;
  a_s = k * dx / dy;
  a_w = k * dy / dx;

  // Boundary conditions
  // Top dirichlet
  a_n.setRow(Ny - 1, 2 * k * dx / dy);
  // Right Neumann
  a_e.setColumn(Nx - 1, 0);
  // Bottom dirichlet
  a_s.setRow(0, 2 * k * dx / dy);
  // Left dirichlet
  a_w.setColumn(0, 2 * k * dy / dx);

  // Center point
  for (unsigned int i = 0; i < Nx; ++i)
    for (unsigned int j = 0; j < Ny; ++j)
      a_p(i, j) = a_n(i, j) + a_e(i, j) + a_s(i, j) + a_w(i, j);
}

void Conduction2D::precomputeRow(unsigned int j) {
  TriDiagonal &A = pre_A_x[j];
  std::vector<double> &b = pre_b_x[j];
  b.resize(Nx);

  // First treat all as an internal volume
  A.addTopRow(a_p(0, j) * w_inv, -a_e(0, j));
  A.addBottomRow(-a_w(Nx - 1, j), a_p(Nx - 1, j) * w_inv);
  for (unsigned int i = 1; i < Nx - 1; ++i)
    A.addMiddleRow(i, -a_w(i, j), a_p(i, j) * w_inv, -a_e(i, j));

  // Left Dirichlet
  b[0] += T_L * a_w(0, j);
  // Top dirichlet
  if (j == Ny - 1)
    for (unsigned int i = 0; i < Nx; ++i)
      b[i] += T_T * a_n(i, j);
  // Bottom dirichlet
  else if (j == 0)
    for (unsigned int i = 0; i < Nx; ++i)
      b[i] += T_B * a_s(i, j);
}

void Conduction2D::precomputeColumn(unsigned int i) {
  TriDiagonal &A = pre_A_y[i];
  std::vector<double> &b = pre_b_y[i];
  b.resize(Ny);

  // First treat all as an internal volume
  A.addTopRow(a_p(i, 0) * w_inv, -a_n(i, 0));
  A.addBottomRow(-a_s(i, Ny - 1), a_p(i, Ny - 1) * w_inv);
  for (unsigned int j = 1; j < Ny - 1; ++j)
    A.addMiddleRow(j, -a_s(i, j), a_p(i, j) * w_inv, -a_n(i, j));

  // Left Dirichlet
  if (i == 0)
    for (unsigned int j = 0; j < Ny; ++j)
      b[j] += T_L * a_w(i, j);
  // Top Dirichlet
  b[Ny - 1] += T_T * a_n(i, Ny - 1);
  // Bottom Dirichlet
  b[0] += T_B * a_s(i, 0);
}

void Conduction2D::precompute() {
  precomputeProperties();
  for (unsigned int i = 0; i < Nx; ++i)
    precomputeColumn(i);
  for (unsigned int j = 0; j < Ny; ++j)
    precomputeRow(j);
}

void Conduction2D::solveRow(unsigned int j) {
  // Copy pre-filled Ax = b for this row
  A_x.copyFrom(pre_A_x[j]);
  b_x.assign(pre_b_x[j].begin(), pre_b_x[j].end());

  // RHS contribution from volumes above and below
  if (j > 0)
    for (unsigned int i = 0; i < Nx; ++i)
      b_x[i] += T(i, j - 1) * a_s(i, j);
  if (j < Ny - 1)
    for (unsigned int i = 0; i < Nx; ++i)
      b_x[i] += T(i, j + 1) * a_n(i, j);

  // Relax, solve, and store solution (which is in b_x)
  for (unsigned int i = 0; i < Nx; ++i)
    b_x[i] += (w_inv - 1.0) * a_p(i, j) * T(i, j);

  A_x.solveTDMA(b_x);
  T.setRow(j, b_x);
}

void Conduction2D::solveColumn(unsigned int i) {
  // Copy pre-filled Ax = b for this row
  b_y.assign(pre_b_y[i].begin(), pre_b_y[i].end());
  A_y.copyFrom(pre_A_y[i]);

  // RHS contribution from volumes left and right
  if (i > 0)
    for (unsigned int j = 0; j < Ny; ++j)
      b_y[j] += T(i - 1, j) * a_w(i, j);
  if (i < Nx - 1)
    for (unsigned int j = 0; j < Ny; ++j)
      b_y[j] += T(i + 1, j) * a_e(i, j);

  // Relax, solve, and store solution (which is in b_y)
  for (unsigned int j = 0; j < Ny; ++j)
    b_y[j] += (w_inv - 1.0) * a_p(i, j) * T(i, j);

  A_y.solveTDMA(b_y);
  T.setColumn(i, b_y);
}

void Conduction2D::sweep() {
  // South to north
  for (int j = 0; j < Ny; ++j)
    solveRow(j);
  // West to east
  for (int i = 0; i < Nx; ++i)
    solveColumn(i);
  // North to south
  for (int j = Ny - 1; j >= 0; --j)
    solveRow(j);
  // East to west
  for (int i = Nx - 1; i >= 0; --i)
    solveColumn(i);
}

double Conduction2D::computeResidual() {
  double R = 0.0, val = 0.0;
  for (unsigned int i = 0; i < Nx; ++i)
    for (unsigned int j = 0; j < Ny; ++j) {
      val = a_p(i, j) * T(i, j) - pre_b_y[i][j];
      if (i > 0)
        val -= a_w(i, j) * T(i - 1, j);
      if (i < Nx - 1)
        val -= a_e(i, j) * T(i + 1, j);
      if (j > 0)
        val -= a_s(i, j) * T(i, j - 1);
      if (j < Ny - 1)
        val -= a_n(i, j) * T(i, j + 1);
      R += std::abs(val);
    }
  return R;
}
