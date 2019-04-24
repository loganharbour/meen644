#include "Problem.h"

using namespace Flow2D;
using namespace std;

void
run(const double Re, const unsigned int Nx, const unsigned int Ny)
{
  const double Lx = 0.6;
  const double Ly = 0.02;
  const double cp = 4183;
  const double k = 0.609;
  const double rho = 998.3;
  const double mu = 0.001002;
  const double q_val = -64;
  const double Sx = 0.06;
  const double Sy = 0.5 * Ly;
  const double T_max = 1.5;
  const double T_w = 0;
  const double u_max = 2 * mu * Re / (rho * Ly);

  // Function to fill a material with v outside of the step and vs inside the step
  auto step_mat = [Sx, Sy](const double v, const double vs, const vector<double> p) -> double {
    // Inside of step
    if (p[0] <= Sx && p[1] >= Sy)
      return vs;
    // Outside of step
    else
      return v;
  };

  // Function to fill u initial condition
  auto u_ic = [u_max, Ly, Sy, Ny](const vector<double> p) -> double {
    // Left side and from 0 < y < Ly / 2
    if ((p[0] == 0) && (p[1] < Sy))
      return u_max * (4 * p[1] / Ly) * (2 - (4 * p[1]) / Ly);
    // Zero otherwise
    else
      return 0;
  };

  // Function to fill v initial condition (all zero)
  auto v_ic = [](const vector<double> p) -> double { return 0; };

  // Function to fill T initial condition
  auto T_ic = [T_max, T_w, Ly, Sy, Ny](const vector<double> p) -> double {
    // Left side and from 0 < y < Ly / 2
    if ((p[0] == 0) && (p[1] < Sy))
      return (T_max - T_w) * (4 * p[1] / Ly) * (2 - (4 * p[1]) / Ly) + T_w;
    // Zero otherwise
    else
      return 0;
  };

  // Function to fill heat flux
  auto q = [q_val, Ly, Sx](const vector<double> p) -> double {
    // Bottom plate or top plate right of step
    if (p[1] == 0 || ((p[1] == Ly) && (p[0] > Sx)))
      return q_val;
    // Zero otherwise
    else
      return 0;
  };

  // Standard inputs
  InputArguments input;
  input.Lx = Lx;
  input.Ly = Ly;
  input.k = [step_mat, k](const vector<double> p) { return step_mat(k, 1E-99, p); };
  input.mu = [step_mat, mu](const vector<double> p) { return step_mat(mu, 1E99, p); };
  input.rho = rho;
  input.cp = cp;
  input.L_ref = Lx;
  input.q = q;
  input.u_ic = u_ic;
  input.v_ic = v_ic;
  input.T_ic = T_ic;
  input.u_ref = u_max;

  Problem problem(Nx, Ny, input);
  problem.run();

  const string prefix = "Re" + to_string((int)Re) + "_Nx" + to_string(Nx) + "_Ny" + to_string(Ny);
  problem.save(Variables::u, "results/" + prefix + "_u.csv");
  problem.save(Variables::v, "results/" + prefix + "_v.csv");
  problem.save(Variables::T, "results/" + prefix + "_T.csv");
}

int
main()
{
  // Problem 1
  for (const unsigned int Ny : {30, 50, 70, 90})
  {
    cout << "Problem 1: Re = 200, 160x" << to_string(Ny) << endl << "  ";
    run(200, 160, Ny);
  }

  // Problem 2
  for (const double Re : {100, 300, 400})
  {
    cout << "Problem 2 - Re = " << to_string((int)Re) << ", 160x60" << endl << "  ";
    run(Re, 160, 70);
  }
}
