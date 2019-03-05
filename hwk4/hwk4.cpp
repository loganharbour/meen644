#include "Flow2D.h"
#include <boost/format.hpp>
#include <map>
#include <sstream>

int
main()
{
  double Re = 1000;
  double Lx = 0.1;
  double Ly = 0.1;
  double mu = 0.001002;
  double rho = 998.3;
  double bc_val = Re * mu / (rho * Lx);
  BoundaryCondition u_BC(bc_val, 0, bc_val, 0);
  BoundaryCondition v_BC(0, 0, 0, 0);

  Flow2D problem(5, 5, Lx, Ly, u_BC, v_BC, rho, mu, false);
  problem.solve();
}
