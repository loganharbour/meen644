#include "TriDiagonal.h"

#include <fstream>
#include <iomanip>
#include <math.h>
#include <string>

void saveCsv(std::vector<double> x, std::string filename)
{
	std::ofstream f;
	f.open(filename);
	for (unsigned int i = 0; i < x.size(); ++i)
		f << std::setprecision(12) << x[i] << std::endl;
	f.close();
}

void solveRod(unsigned int N)
{
	// Initialize geometry and material properties
	double L = 1;       // [m]
	double dx = L / N;  // [m]
	double r = 0.025;   // [m]
	double k = 400.0;   // [W/m-C]
	double h = 0.5;     // [W/m^2-C]
	double T0 = 100.0;  // [C]
	double Tinf = 25.0; // [C]

	// Initialize system A theta = b
	TriDiagonal A(N - 1);
	std::vector<double> theta(N - 1);
	std::vector<double> b(N - 1);

	// Constant coefficients
	double a_p = 2.0 * h * dx / (k * r) + 2.0 / dx;
	double a_w = 1.0 / dx;
	double a_e = 1.0 / dx;

	// Fill system
	A.setTopRow(a_p, -a_e);
	b[0] = a_w * (T0 - Tinf);
	for (unsigned int i = 1; i < N - 2; ++i)
		A.setMiddleRow(i, -a_w, a_p, -a_e);
	A.setBottomRow(-1, 1 + h * dx / (2 * k));

	// Solve system in place
	A.solveTDMA(b, theta);

	// Add in Dirichlet left node and save to csv
	std::vector<double> T(N);
	T[0] = T0;
	for (unsigned int i = 0; i < theta.size(); ++i)
		T[i + 1] = theta[i] + Tinf;
	saveCsv(T, "../results/result_" + std::to_string(N) + ".csv");
}

int main()
{
	solveRod(6);
	solveRod(11);
	solveRod(21);
	solveRod(41);
	solveRod(81);

	return 0;
}
