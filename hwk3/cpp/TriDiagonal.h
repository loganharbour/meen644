#ifndef TRIDIAGONAL_H
#define TRIDIAGONAL_H

#include "Base.h"

/**
 * Class that holds a tri-diagonal matrix and is able to perform TDMA in place
 * with a given RHS.
 */
class TriDiagonal {
public:
  TriDiagonal(unsigned int N, double val = 0)
      : N(N), A(N, val), B(N, val), C(N - 1, val) {}

  // Gets the value of the (i, j) entry
  const double operator()(unsigned int i, unsigned int j) {
    assert(i < N && j > i - 2 && j < i + 2);
    if (j == i - 1)
      return A[i];
    else if (j == i)
      return B[i];
    else if (j == i + 1)
      return C[i];
    else {
      std::cerr << "( " << i << "," << j << ") out of TriDiagonal system";
      std::terminate();
    }
  }

  // Adders for the top, middle, and bottom rows
  void addTopRow(double b, double c) {
    B[0] += b;
    C[0] += c;
  }
  void addMiddleRow(unsigned int i, double a, double b, double c) {
    assert(i < N - 1 && i != 0);
    A[i] += a;
    B[i] += b;
    C[i] += c;
  }
  void addBottomRow(double a, double b) {
    A[N - 1] += a;
    B[N - 1] += b;
  }

  // Adders for the left, main, and right diagonals
  void addLeft(unsigned int i, double val) {
    assert(i < N && i != 0);
    A[i] += val;
  }
  void addMain(unsigned int i, double val) {
    assert(i < N);
    B[i] += val;
  }
  void addRight(unsigned int i, double val) {
    assert(i < N - 1);
    C[i] += val;
  }

  void copyFrom(TriDiagonal &from) {
    assert(from.getN() == N);
    A.assign(from.getA().begin(), from.getA().end());
    B.assign(from.getB().begin(), from.getB().end());
    C.assign(from.getC().begin(), from.getC().end());
  }

  const std::vector<double> &getA() { return A; }
  const std::vector<double> &getB() { return B; }
  const std::vector<double> &getC() { return C; }
  unsigned int getN() { return N; }

  // Solves the system Ax = d in place where d eventually stores the solution
  void solveTDMA(std::vector<double> &d) {
    // Forward sweep
    double tmp = 0;
    for (unsigned int i = 1; i < N; ++i) {
      tmp = A[i] / B[i - 1];
      B[i] -= tmp * C[i - 1];
      d[i] -= tmp * d[i - 1];
    }

    // Backward sweep
    d[N - 1] /= B[N - 1];
    for (unsigned int i = N - 2; i != std::numeric_limits<unsigned int>::max();
         --i) {
      d[i] -= C[i] * d[i + 1];
      d[i] /= B[i];
    }
  }

protected:
  // Matrix size (N x N)
  unsigned int N;
  // Left/main/right diagonal storage
  std::vector<double> A, B, C;
};

#endif /* TRIDIAGONAL_H */
