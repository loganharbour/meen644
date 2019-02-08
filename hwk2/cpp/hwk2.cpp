#include "Base.h"
#include "TriDiagonal.h"

/**
 * Solves the 1D heat-conduction (with convection) problem with N nodes.
 */
void solveRod(unsigned int N) {
  // Initialize geometry, material properties, and constant coefficients
  double L = 1,         // [m]
      dx = L / (N - 1), // [m]
      d = 0.05,         // [m]
      k = 400.0,        // [W/m-C]
      h = 0.5,          // [W/m^2-C]
      T0 = 100.0,       // [C]
      Tinf = 25.0,      // [C]
      m = 4.0 * h / (k * d),
      a_p = m * dx + 2.0 / dx,
      a_w = 1.0 / dx,
      a_e = 1.0 / dx;

  // Initialize system A theta = b
  TriDiagonal A(N - 1);
  std::vector<double> theta(N - 1), b(N - 1);

  // Fill system
  A.addTopRow(a_p, -a_e);
  b[0] = a_w * (T0 - Tinf);
  for (unsigned int i = 1; i < N - 2; ++i)
    A.addMiddleRow(i, -a_w, a_p, -a_e);
  A.addBottomRow(-a_w, m * dx / 2.0 + h / k + a_e);

  // Solve system in place and save
  A.solveTDMA(b, theta);
  saveVectorCsv(theta, "../results/theta_" + std::to_string(N) + ".csv");
}

int main() {
  // Run each requested case
  for (unsigned int N : {6, 11, 21, 41, 81})
    solveRod(N);

  return 0;
}
