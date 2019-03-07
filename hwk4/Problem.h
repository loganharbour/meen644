#ifndef PROBLEM_H
#define PROBLEM_H

#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>

#include "Variable.h"

namespace Flow2D
{

using namespace std;

struct InputArguments
{
  double Lx, Ly;
  BoundaryCondition u_bc, v_bc;
  double L_ref, u_ref;
  double rho, mu;
  bool debug = false;
  double alpha_p = 0.7;
  double alpha_uv = 0.5;
  unsigned int max_its = 100000;
  double tol = 1.0e-6;
};

class Problem
{
public:
  Problem(const unsigned int Nx, const unsigned int Ny, const InputArguments & input);

  void run();

  // Public access to printing and saving variable results
  void print(const Variables var,
             const string prefix = "",
             const bool newline = false,
             const unsigned int pr = 5) const
  {
    variables.at(var).print(prefix, newline, pr);
  }
  void save(const Variables var, const string filename) const { variables.at(var).save(filename); }

private:
  // Problem_corrections.cpp
  void correct();
  void pCorrect();
  void uCorrect();
  void vCorrect();

  // Problem_coefficients.cpp
  void fillCoefficients(Variable & var);
  void pcCoefficients();
  void uCoefficients();
  void vCoefficients();
  void velocityCoefficients(Coefficients & a,
                            const Coefficients & D,
                            const Coefficients & F,
                            const double & b);

  // Problem_residuals.cpp
  void computeResiduals();
  double pResidual() const;
  double velocityResidual(const Variable & var) const;

  // Problem_solvers.cpp
  void solve();
  void solve(Variable & var);
  void sweepColumns(Variable & var, const bool west_east = true);
  void sweepRows(Variable & var, const bool south_north = true);
  void sweepColumn(const unsigned int i, Variable & var);
  void sweepRow(const unsigned int j, Variable & var);
  void solveVelocities();

protected:
  // Number of pressure CVs
  const unsigned int Nx, Ny;

  // Geometry [m]
  const double Lx, Ly, dx, dy;
  // Material properties
  const double rho, mu;
  // Residual references
  const double L_ref, u_ref;

  // Enable debug mode (printing extra output)
  const bool debug;

  // Maximum iterations
  const unsigned int max_its;
  // Iteration tolerance
  const double tol;
  // Pressure relaxation
  const double alpha_p;
  // Number of iterations completed
  unsigned int iterations = 0;

  // Variables
  Variable u, v, pc, p;
  // Variable map
  map<const Variables, const Variable &> variables;

  // Whether or not we converged
  bool converged = false;
};

} // namespace Flow2D
#endif /* PROBLEM_H */
