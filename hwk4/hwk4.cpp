#include "Problem.h"

using namespace Flow2D;

int
main()
{
  // Problem wide constants
  double L = 0.1;
  double Re = 400;
  double rho = 998.3;
  double mu = 1.002e-3;
  double bc_val = Re * mu / (rho * 0.1);

  // Standard inputs
  InputArguments input;
  input.Lx = L;
  input.Ly = L;
  input.mu = mu;
  input.rho = rho;
  input.u_ref = bc_val;
  input.L_ref = L;
  input.u_bc = BoundaryCondition(0, bc_val, 0, bc_val);
  input.v_bc = BoundaryCondition(0, 0, 0, 0);

  // Problem 2: change to top plate BC
  input.u_bc = BoundaryCondition(bc_val, 0, 0, 0);
  input.v_bc = BoundaryCondition(0, 0, 0, 0);
  std::cout << "Problem 2, top plate BC" << std::endl;
  for (unsigned int N : {8, 16, 32, 64, 128})
  {
    std::cout << "  N = " << N << "x" << N << " - ";
    Problem problem(N, N, input);
    problem.run();
    problem.save(Variables::u, "results/p2/" + to_string(N) + "_u.csv");
    problem.save(Variables::v, "results/p2/" + to_string(N) + "_v.csv");
  }
}
