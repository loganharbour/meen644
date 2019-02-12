#ifndef TRIDIAGONAL_H
#define TRIDIAGONAL_H

#include "Base.h"

/**
 * Class that holds a tri-diagonal matrix and is able to perform TDMA in place
 * with a given RHS and solution vector.
 */
class TriDiagonal {
public:
  TriDiagonal(unsigned int N) : N(N), A(N), B(N), C(N - 1) {}

  // Adds val to (i, j)
  void add(unsigned int i, unsigned int j, double val) {
    assert(i < N && j > i - 2 && j < i + 2);
    if (j == i - 1)
      A[i] += val;
    else if (j == i)
      B[i] += val;
    else if (j == i + 1)
      C[i] += val;
    else {
      std::cerr << "( " << i << "," << j << ") out of TriDiagonal system";
      std::terminate();
    }
  }

  // Adders for the top, interior, and bottom rows
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

  // Reset the matrix to zero
  void clear() {
    clearVector(A);
    clearVector(B);
    clearVector(C);
  }

  // Saves the three diagonal vectors in separate csv files
  void save(const std::string file_prefix) {
    saveVectorCsv(A, file_prefix + "A.csv");
    saveVectorCsv(B, file_prefix + "B.csv");
    saveVectorCsv(C, file_prefix + "C.csv");
  }

  // Solves the system Ax = d in place where A is the matrix held by this class
  void solveTDMA(std::vector<double> &d, std::vector<double> &x) {
    double w = 0;

    // Forward sweep
    for (unsigned int i = 1; i < N; ++i) {
      w = A[i] / B[i - 1];
      B[i] -= w * C[i - 1];
      d[i] -= w * d[i - 1];
    }

    // Backward substitution
    x[N - 1] = d[N - 1] / B[N - 1];
    for (unsigned int i = N - 2; i != std::numeric_limits<unsigned int>::max();
         --i)
      x[i] = (d[i] - C[i] * x[i + 1]) / B[i];
  }

  // Getters for the diagonal vectors
  const std::vector<double> &a() { return A; }
  const std::vector<double> &b() { return B; }
  const std::vector<double> &c() { return C; }

protected:
  // Matrix size (N x N)
  unsigned int N;
  // Left/main/right diagonal storage
  std::vector<double> A, B, C;
};

#endif /* TRIDIAGONAL_H */
