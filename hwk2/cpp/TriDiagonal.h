#include "Base.h"

/**
 * Class that holds a tri-diagonal matrix and is able to perform TDMA in place
 * with a given RHS and solution vector.
 */
class TriDiagonal {
public:
  TriDiagonal(unsigned int N) : _N(N), _a(N), _b(N), _c(N - 1) {}

  // Adds val to (i, j)
  void add(unsigned int i, unsigned int j, double val) {
    assert(i < _N && j > i - 2 && j < i + 2);
    if (j == i - 1)
      _a[i] += val;
    else if (j == i)
      _b[i] += val;
    else if (j == i + 1)
      _c[i] += val;
    else {
      std::cerr << "( " << i << "," << j << ") out of TriDiagonal system";
      std::terminate();
    }
  }

  // Adders for the top, interior, and bottom rows
  void addTopRow(double b, double c) {
    _b[0] += b;
    _c[0] += c;
  }
  void addMiddleRow(unsigned int i, double a, double b, double c) {
    assert(i < _N - 1 && i != 0);
    _a[i] += a;
    _b[i] += b;
    _c[i] += c;
  }
  void addBottomRow(double a, double b) {
    _a[_N - 1] += a;
    _b[_N - 1] += b;
  }

  // Adders for the left, main, and right diagonals
  void addLeft(unsigned int i, double val) {
    assert(i < _N && i != 0);
    _a[i] += val;
  }
  void addMain(unsigned int i, double val) {
    assert(i < _N);
    _b[i] += val;
  }
  void addRight(unsigned int i, double val) {
    assert(i < _N - 1);
    _c[i] += val;
  }

  // Reset the matrix to zero
  void reset() {
    std::fill(_a.begin(), _a.end(), 0);
    std::fill(_b.begin(), _b.end(), 0);
    std::fill(_c.begin(), _c.end(), 0);
  }

  // Saves the three diagonal vectors in separate csv files
  void saveCsv(const std::string file_prefix) {
    saveVectorCsv(_a, file_prefix + "_a.csv");
    saveVectorCsv(_b, file_prefix + "_b.csv");
    saveVectorCsv(_c, file_prefix + "_c.csv");
  }

  // Solves the system Ax = d in place where A is the matrix held by this class
  void solveTDMA(std::vector<double> &d, std::vector<double> &x) {
    double w = 0;

    // Forward sweep
    for (unsigned int i = 1; i < _N; ++i) {
      w = _a[i] / _b[i - 1];
      _b[i] -= w * _c[i - 1];
      d[i] -= w * d[i - 1];
    }

    // Backward substitution
    x[_N - 1] = d[_N - 1] / _b[_N - 1];
    for (unsigned int i = _N - 2; i != std::numeric_limits<unsigned int>::max();
         --i)
      x[i] = (d[i] - _c[i] * x[i + 1]) / _b[i];
  }

  // Getters for the diagonal vectors
  const std::vector<double> &a() { return _a; }
  const std::vector<double> &b() { return _b; }
  const std::vector<double> &c() { return _c; }

protected:
  // Matrix size (_N x _N)
  unsigned int _N;
  // Left/main/right diagonal storage
  std::vector<double> _a, _b, _c;
};
