#include <vector>

class TriDiagonal
{
public:
	TriDiagonal(unsigned int N) : N(N),	A(N), B(N), C(N)
	{
	}

	void setTopRow(double b, double c)
	{
		B[0] = b;
		C[0] = c;
	}

	void setMiddleRow(unsigned int i, double a, double b, double c)
	{
		A[i] = a;
		B[i] = b;
		C[i] = c;
	}

	void setBottomRow(double a, double b)
	{
		A[N - 1] = a;
		B[N - 1] = b;
	}

	void solveTDMA(std::vector<double> & D, std::vector<double> & X)
	{
		double w = 0;

		// Forward sweep
		for (unsigned int i = 1; i < N; ++i)
		{
			w = A[i] / B[i - 1];
			B[i] -= w * C[i - 1];
			D[i] -= w * D[i - 1];
		}

		// Backward sweep
		X[N - 1] = D[N - 1] / B[N - 1];
		for (unsigned int i = N - 2; i != std::numeric_limits<unsigned int>::max(); --i)
			X[i] = (D[i] - C[i] * X[i + 1]) / B[i];
	}


protected:
	// Matrix size (N x N)
	unsigned int N;

	// The three diagonals
	std::vector<double> A;
	std::vector<double> B;
	std::vector<double> C;
};
