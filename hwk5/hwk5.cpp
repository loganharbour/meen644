#include "Problem.h"

using namespace Flow2D;

int
main()
{
  // Problem wide constants
  double Lx = 2;
  double Ly = 0.02;
  double cp = 4183;
  double k = 0.609;
  double rho = 998.3;
  double mu = 0.001002;
  double Re = 100;
  double u_bc_val = Re * mu / (2 * rho * Ly);
  double T_bc_val = 25;
  double q_bc_val = 500;

  // Standard inputs
  InputArguments input;
  input.Lx = Lx;
  input.Ly = Ly;
  input.cp = cp;
  input.k = k;
  input.mu = mu;
  input.rho = rho;
  input.u_ref = u_bc_val;
  input.L_ref = Lx;
  input.T_left = T_bc_val;
  input.u_bc = BoundaryCondition(0, 0, 0, u_bc_val);
  input.v_bc = BoundaryCondition(0, 0, 0, 0);
  input.q_top_bot = q_bc_val;

  // Problem 3: check symmetry
  std::cout << "Problem 3: check symmetry" << std::endl;
  Problem problem3(10, 5, input);
  problem3.run();
  problem3.print(Variables::u, "u =");
  problem3.print(Variables::p, "p =");
  problem3.print(Variables::T, "T = ", true);

  // Problem 4: 180x54 grid
  std::cout << "Problem 4: 180 x 54 grid" << std::endl;
  Problem problem4(180, 54, input);
  problem4.run();
  problem4.save(Variables::u, "results/u.csv");
  problem4.save(Variables::v, "results/v.csv");
  problem4.save(Variables::T, "results/T.csv");
}
