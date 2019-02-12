#include "Base.h"
#include "Solution.h"
#include "TriDiagonal.h"

class Conduction2D {
public:
  Conduction2D(unsigned int Nx, unsigned int Ny, double w, double Lx = 0.5,
               double Ly = 0.5, double k = 386.0, double T_L = 50.0,
               double T_T = 100.0, double T_B = 50.0, double tol = 1e-5)
      : // Interior nodal points
        Nx(Nx), Ny(Ny),
        // Domain size
        Lx(Lx), Ly(Ly), dx(Lx / Nx), dy(Ly / Ny), // [m]
        // Material properties
        k(k), // [W / m - K]
        a_n(k * dx / dy), a_e(k * dy / dx), a_s(k * dx / dy), a_w(k * dy / dx),
        a_p(a_n + a_e + a_s + a_w),
        // Boundary conditions
        T_L(T_L), T_T(T_T), T_B(T_B),
        // Relaxation coefficient
        w(w),
        // Initialize arrays
        T(Nx, Ny), Tx(Nx), Ty(Ny), bx(Nx), by(Ny), temp_x(Nx), temp_y(Nx),
        Ax(Nx), Ay(Ny) {}

  void solveRow(unsigned int row) {
    assert(row < Ny);
    Ax.clear();
    clearVector(bx);

    // Treat every row the same on the LHS and correct later
    Ax.addTopRow(a_p, -a_e);
    Ax.addBottomRow(-a_w, a_p);
    for (unsigned int i = 1; i < Nx - 1; ++i)
      Ax.addMiddleRow(i, -a_w, a_p, -a_e);

    // RHS contribution from lagged above
    if (row > 0) {
      T.getRow(row - 1, temp_x);
      for (unsigned int i = 0; i < Nx; ++i)
        bx[i] += temp_x[i] * a_s;
    }
    // RHS contribution from newly updated below
    if (row < Ny - 1) {
      T.getRow(row + 1, temp_x);
      for (unsigned int i = 0; i < Nx; ++i)
        bx[i] += temp_x[i] * a_n;
    }

    // Left Dirichlet
    bx[0] += 2 * T_L * a_w;
    Ax.addMain(0, a_w);
    // Top Dirichlet
    if (row == Ny - 1)
      for (unsigned int i = 0; i < Nx; ++i) {
        bx[i] += 2 * T_T * a_n;
        Ax.addMain(i, a_n);
      }
    // Bottom Dirichlet
    else if (row == 0)
      for (unsigned int i = 0; i < Nx; ++i) {
        bx[i] += 2 * T_B * a_s;
        Ax.addMain(i, a_s);
      }
    // Right Neumann
    Ax.addMain(Nx - 1, -a_e);

    // Solve and store
    Ax.solveTDMA(bx, Tx);
    T.setRow(row, Tx);
  }

  void solveColumn(unsigned int col) {
    assert(col < Nx);
    Ay.clear();
    clearVector(by);

    // Treat every column the same on the LHS and correct later
    Ay.addTopRow(a_p, -a_n);
    Ay.addBottomRow(-a_s, a_p);
    for (unsigned int j = 1; j < Ny - 1; ++j)
      Ay.addMiddleRow(j, -a_s, a_p, -a_n);

    // RHS contribution from newly updated left
    if (col > 0) {
      T.getColumn(col - 1, temp_y);
      for (unsigned int j = 0; j < Ny; ++j)
        by[j] += temp_y[j] * a_w;
    }
    // RHS contribution from lagged right
    if (col < Nx - 1) {
      T.getColumn(col + 1, temp_y);
      for (unsigned int j = 0; j < Ny; ++j)
        by[j] += temp_y[j] * a_n;
    }

    // Left Dirichlet
    if (col == 0)
      for (unsigned int j = 0; j < Ny; ++j) {
        by[j] += 2 * T_L * a_w;
        Ay.addMain(j, a_w);
      }
    // Top Dirichlet
    by[Ny - 1] += 2 * T_T * a_n;
    Ay.addMain(Ny - 1, a_n);
    // Bottom Dirichlet
    by[0] += 2 * T_B * a_s;
    Ay.addMain(0, a_s);
    // Right Neumann
    if (col == Ny - 1)
      for (unsigned int j = 0; j < Ny; ++j)
        Ay.addMain(j, -a_e);

    // Solve and store
    Ay.solveTDMA(by, Ty);
    T.setColumn(col, Ty);
  }

  void sweep() {
    for (unsigned int i = 0; i < Nx; ++i)
      solveRow(i);
    for (unsigned int j = 0; j < Ny; ++j)
      solveColumn(j);
    std::cout << computeResidual() << std::endl;
  }

  const double computeResidual() {
    double R = 0, val = 0;
    for (unsigned int i = 0; i < Nx; ++i)
      for (unsigned int j = 0; j < Ny; ++j) {
        val = a_p * T.get(i, j);
        if (j < Ny - 1)
          val -= a_n * T.get(i, j + 1);
        else
          val -= a_n * T_T;
        if (i < Nx - 1)
          val -= a_e * T.get(i + 1, j);
        if (j > 0)
          val -= a_s * T.get(i, j - 1);
        else
          val -= a_s * T_B;
        if (i > 0)
          val -= a_w * T.get(i - 1, j);
        else
          val -= a_w * T_L;
        R += std::abs(val);
      }
    return R;
  }

  void run() {
    for (unsigned int l = 0; l < 100; ++l)
      sweep();
  }

  const void save(std::string filename) { T.save(filename); }

protected:
  // Number of interior nodal points in the x and y-dimensions
  const unsigned int Nx, Ny;

  // Material properties, geometric properties, boundary conditions
  const double Lx, Ly, dx, dy, k;
  const double a_n, a_e, a_s, a_w, a_p;
  const double T_L, T_T, T_B;

  // Relaxation coefficient
  const double w;

  // Current temperature solution
  Solution T;

  // Arrays for the iterative solves
  std::vector<double> Tx, Ty, bx, by, temp_x, temp_y;
  TriDiagonal Ax, Ay;
};
