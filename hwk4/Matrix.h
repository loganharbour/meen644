#ifndef MATRIX
#define MATRIX

// #define NDEBUG
#include <cassert>
#include <vector>

/**
 * Class that holds a N x M matrix with common matrix operations.
 */
template <typename T>
class Matrix {
public:
  Matrix(unsigned int N, unsigned int M)
      : N(N), M(M), A(N, std::vector<T>(M)) {}

  // Const operator for getting the (i, j) element
  const T &operator()(unsigned int i, unsigned int j) const {
    assert(i < N && j < M);
    return A[i][j];
  }
  // Operator for getting the (i, j) element
  T &operator()(unsigned int i, unsigned int j) {
    assert(i < N && j < M);
    return A[i][j];
  }
  // Operator for setting the entire matrix to a value
  void operator=(T v) {
    for (unsigned int j = 0; j < M; ++j)
      setRow(j, v);
  }

  // Saves the matrix in csv format
  void save(const std::string filename, unsigned int precision = 12) const {
    std::ofstream f;
    f.open(filename);
    for (unsigned int j = 0; j < M; ++j) {
      for (unsigned int i = 0; i < N; ++i) {
        if (i > 0)
          f << ",";
        f << std::setprecision(precision) << A[i][j];
      }
      f << std::endl;
    }
    f.close();
  }

  // Set the j-th row to v
  void setRow(unsigned int j, T v) {
    assert(j < M);
    for (unsigned int i = 0; i < N; ++i)
      A[i][j] = v;
  }
  // Set the i-th column to v
  void setColumn(unsigned int i, T v) {
    assert(i < N);
    for (unsigned int j = 0; j < M; ++j)
      A[i][j] = v;
  }

  // Set the j-th row to vector v
  void setRow(unsigned int j, std::vector<T> &v) {
    assert(j < M && v.size() == N);
    for (unsigned int i = 0; i < N; ++i)
      A[i][j] = v[i];
  }
  // Set the i-th column to vector v
  void setColumn(unsigned int i, std::vector<T> &v) {
    assert(i < N && v.size() == M);
    for (unsigned int j = 0; j < M; ++j)
      A[i][j] = v[j];
  }

private:
  // The size of this matrix
  const unsigned int N, M;

  // Matrix storage
  std::vector<std::vector<T> > A;
};

#endif /* MATRIX_H */
