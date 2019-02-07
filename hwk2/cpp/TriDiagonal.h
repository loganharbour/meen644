#include <vector>

/**
 * Class that holds a tri-diagonal matrix and is able to perform TDMA with a
 * given RHS vector.
 */
class TriDiagonal
{
public:
	TriDiagonal(unsigned int N) : _N(N), _a(N), _b(N), _c(N - 1)
	{
	}

	// Sets the top row diagonal to b and off diagonal to c
	void setTopRow(double b, double c)
	{
		_b[0] = b;
		_c[0] = c;
	}

	// Sets an internal row's left off-diagonal to a, its diagonal to b, and its
	// right off-diagonal to c
	void setMiddleRow(unsigned int i, double a, double b, double c)
	{
		_a[i] = a;
		_b[i] = b;
		_c[i] = c;
	}

	// Sets the bottom row off diagonal to a and diagonal to b
	void setBottomRow(double a, double b)
	{
		_a[_N - 1] = a;
		_b[_N - 1] = b;
	}

	// Solves the system Ax = d in place where A is the matrix held by this class
	void solveTDMA(std::vector<double> & d, std::vector<double> & x)
	{
		double w = 0;

		// Forward sweep
		for (unsigned int i = 1; i < _N; ++i)
		{
			w = _a[i] / _b[i - 1];
			_b[i] -= w * _c[i - 1];
			d[i] -= w * d[i - 1];
		}

		// Backward substitution
		x[_N - 1] = d[_N - 1] / _b[_N - 1];
		for (unsigned int i = _N - 2; i != std::numeric_limits<unsigned int>::max(); --i)
			x[i] = (d[i] - _c[i] * x[i + 1]) / _b[i];
	}

	// Getters for the diagonal vectors
	const std::vector<double> & a() { return _a; }
	const std::vector<double> & b() { return _b; }
	const std::vector<double> & c() { return _c; }

protected:
	// Matrix size (_N x _N)
	unsigned int _N;
	// Left off-diagonal storage
	std::vector<double> _a;
	// Diagonal storage
	std::vector<double> _b;
	// Right off-diagonal storage
	std::vector<double> _c;
};
