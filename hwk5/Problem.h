#ifndef PROBLEM_H
#define PROBLEM_H

#include <cmath>
#include <ctime>
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
  double T_left, q_top_bot;
  double L_ref, u_ref;
  double cp, k, mu, rho;
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
  void correctMain();
  void correctAux();
  void pCorrect();
  void pBCCorrect();
  void TBCCorrect();
  void uCorrect();
  void uBCCorrect();
  void vCorrect();

  // Problem_coefficients.cpp
  void fillCoefficients(const Variable & var);
  void pcCoefficients();
  void TCoefficients();
  void uCoefficients();
  void vCoefficients();
  void fillPowerLaw(Coefficients & a,
                    const Coefficients & D,
                    const Coefficients & F,
                    const double & b = 0);

  // Problem_residuals.cpp
  void computeMainResiduals();
  void computeAuxResiduals();
  double pResidual() const;
  double TResidual() const;
  double velocityResidual(const Variable & var) const;

  // Problem_solvers.cpp
  void solveMain();
  void solveAux();
  void solve(Variable & var);
  void sweepColumns(Variable & var, const bool west_east = true);
  void sweepRows(Variable & var, const bool south_north = true);
  void sweepColumn(const unsigned int i, Variable & var);
  void sweepRow(const unsigned int j, Variable & var);
  void solveVelocities();

  // Quicker v^5 for velocityCoefficients() (yes, it's actually much faster...)
  static const double pow5(const double & v) { return v * v * v * v * v; }

protected:
  // Number of pressure CVs
  const unsigned int Nx, Ny;

  // Geometry [m]
  const double Lx, Ly, dx, dy;
  // Material properties
  const double cp, k, mu, rho;
  // Residual references
  const double L_ref, u_ref;
  // Other boundary conditions
  const double q_top_bot;
  // Mass inflow
  const double m_in;

  // Enable debug mode (printing extra output)
  const bool debug;

  // Maximum iterations
  const unsigned int max_its;
  // Iteration tolerance
  const double tol;
  // Pressure relaxation
  const double alpha_p;
  // Number of iterations completed
  unsigned int main_iterations = 0;
  unsigned int aux_iterations = 0;

  // Variables
  Variable u, v, pc, p, T;
  // Variable map
  map<const Variables, const Variable &> variables;

  // Whether or not we converged
  bool main_converged = false;
  bool converged = false;
  // Run start time
  clock_t start;
};

} // namespace Flow2D
#endif /* PROBLEM_H */
