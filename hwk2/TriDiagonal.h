#ifndef TRIDIAGONAL_H
#define TRIDIAGONAL_H

#include <iomanip>
#include <fstream>
#include <iostream>
#include <vector>

class TriDiagonal
{
public:
	TriDiagonal(unsigned int N) :
		_N(N),
		_T(_N * 3 - 2)
	{
	}

	void setTopRow(double d, double r)
	{
		_T[0] = d;
		_T[1] = r;
	}

	void setMiddleRow(unsigned int i, double l, double d, double r)
	{
		_T[i * 3 - 1] = l;
		_T[i * 3] = d;
		_T[i * 3 + 1] = r;
	}

	void setBottomRow(double l, double d)
	{
		_T[_N * 3 - 4] = l;
		_T[_N * 3 - 3] = d;
	}

	void solveTDMA(std::vector<double> & b, std::vector<double> & x)
	{
		// Location of the diagonal
		unsigned int d;

		// Forward sweep, first row
		_T[1] /= _T[0];
		b[0] /= _T[0];

		// Forward sweep, middle rows
		for (unsigned int i = 1; i < _N - 1; ++i)
		{
			d = i * 3;
			_T[d + 1] /= _T[d] - _T[d - 1] * _T[d - 2];
			b[i] = (b[i] - _T[d - 1] * b[i - 1]) / (_T[d] - _T[d - 1] * _T[d - 2]);
		}

		// Forward sweep, last row
		d = (_N - 1) * 3;
		b[_N - 1] = (b[_N - 1] - _T[d - 1] * b[_N - 2]) / (_T[d] - _T[d - 1] * _T[d - 2]);

		// Backward sweep, last row
		x[_N - 1] = b[_N - 1];

		// Backward sweep, remaining rows
		for (unsigned int i = _N - 2; i != std::numeric_limits<unsigned int>::max(); --i)
			x[i] = b[i] - _T[i * 3 + 1] * x[i + 1];
	}

protected:
	// Matrix size
	unsigned int _N;

	// Single vector that stores the matrix
	std::vector<double> _T;
};

#endif
