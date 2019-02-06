#include "TriDiagonal.h"
#include <time.h>
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
	double L = 1; // [cm]
	double dx = L / N; // [cm]
	double diameter = 0.05; // [cm]
	double k = 4.0; // [W/cm-C]
	double T0 = 100.0; // [C]
	double Tinf = 25.0; // [C]
	double h = 0.5; // [W/m^2-C]
 	double m = 0.5 * h * diameter;

	TriDiagonal A(N - 1);
	std::vector<double> x(N - 1);
	std::vector<double> T(N);
	std::vector<double> b(N - 1);

	// These coefficients are constant in space
	double a_p = m * dx + 2.0 / dx;
	double a_w = 1.0 / dx;
	double a_e = 1.0 / dx;

	// Left control volume
	A.setTopRow(a_p, -a_e);
	b[0] = a_w * (T0 - Tinf);
	// Interior control volumes
	for (unsigned int i = 1; i < N - 2; ++i)
		A.setMiddleRow(i, -a_w, a_p, -a_e);
	// Right control volume
	A.setBottomRow(-1, 1 + h * dx / (2 * k));

	A.solveTDMA(b, x);

	T[0] = T0;
	for (unsigned int i = 0; i < x.size(); ++i)
		T[i + 1] = x[i] + Tinf;
	saveCsv(T, "result_" + std::to_string(N) + ".csv");
}

int main()
{
	solveRod(6);
	solveRod(11);
	solveRod(21);
	solveRod(41);
	solveRod(81);
	solveRod(10000);
	solveRod(100000);

	return 0;
}