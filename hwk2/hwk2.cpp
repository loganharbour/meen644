#include "TriDiagonal.h"
#include <time.h>

#include <math.h>

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
	double L = 1.0; // [cm]
	double dx = L / N; // [cm]
	double diameter = 0.05; // [cm]
	double k = 400.0;
	double T1 = 100.0; // [C]
	double Tinf = 25.0; // [C]
	double h = 0.5; // [W/m^2-C]
 	double m = 0.5 * h * diameter;

	TriDiagonal M(N - 1);
	std::vector<double> x(N - 1);
	std::vector<double> b(N - 1);

	// These coefficients are constant in space
	double a_p = m * dx + 2.0 / dx;
	double a_w = 1.0 / dx;
	double a_e = 1.0 / dx;

	// Left control volume
	M.setTopRow(a_p, -a_e);
	b[0] = a_w * T1;
	// Interior control volumes
	for (unsigned int i = 1; i < N - 1; ++i)
		M.setMiddleRow(i, -a_w, a_p, -a_e);
	// Right control volume
	M.setBottomRow(-1, 1 - h * dx / (2 * k));

	M.solveTDMA(b, x);
	for (unsigned int i = 0; i < x.size(); ++i)
		x[i] -= Tinf;
	saveCsv(x, "test.csv");

}

int main()
{
	solveRod(100);

	return 0;
}
