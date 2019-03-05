#ifndef TRIDIAGONAL_H
#define TRIDIAGONAL_H

#define NDEBUG
#include <cassert>

/**
 * Class that holds a tri-diagonal matrix and is able to perform TDMA in place
 * with a given RHS.
 */
template <typename T>
class TriDiagonal
{
public:
  TriDiagonal(unsigned int N, T v = 0) : N(N), A(N, v), B(N, v), C(N - 1, v) {}

  // Operator for setting the entire matrix to a value
  void operator=(TriDiagonal & from)
  {
    assert(from.getN() == N);
    A = from.getA();
    B = from.getB();
    C = from.getC();
  }

  // Gets the value of the (i, j) entry
  const T operator()(unsigned int i, unsigned int j) const
  {
    assert(i < N && j > i - 2 && j < i + 2);
    if (j == i - 1)
      return A[i];
    else if (j == i)
      return B[i];
    else if (j == i + 1)
      return C[i];
    else
    {
      std::cerr << "( " << i << "," << j << ") out of TriDiagonal system";
      std::terminate();
    }
  }

  // Adders for the top, middle, and bottom rows
  void addTopRow(T b, T c)
  {
    B[0] += b;
    C[0] += c;
  }
  void addMiddleRow(unsigned int i, T a, T b, T c)
  {
    assert(i < N - 1 && i != 0);
    A[i] += a;
    B[i] += b;
    C[i] += c;
  }
  void addBottomRow(T a, T b)
  {
    A[N - 1] += a;
    B[N - 1] += b;
  }

  // Setters for the top, middle, and bottom rows
  void setTopRow(T b, T c)
  {
    B[0] = b;
    C[0] = c;
  }
  void setMiddleRow(unsigned int i, T a, T b, T c)
  {
    assert(i < N - 1 && i != 0);
    A[i] = a;
    B[i] = b;
    C[i] = c;
  }
  void setBottomRow(T a, T b)
  {
    A[N - 1] = a;
    B[N - 1] = b;
  }

  // Getters for the raw vectors
  const std::vector<T> & getA() const { return A; }
  const std::vector<T> & getB() const { return B; }
  const std::vector<T> & getC() const { return C; }

  // Getter for the size
  unsigned int getN() { return N; }

  // Prints the matrix
  void print(const unsigned int pr = 6)
  {
    for (unsigned int i = 0; i < N; ++i)
      std::cout << std::showpos << std::scientific << std::setprecision(pr) << (i > 0 ? A[i] : 0)
                << " " << B[i] << " " << (i < N - 1 ? C[i] : 0) << std::endl;
  }
  // Saves the matrix in csv format
  void save(const std::string filename, unsigned int precision = 12) const
  {
    std::ofstream f;
    f.open(filename);
    for (unsigned int i = 0; i < N; ++i)
    {
      if (i > 0)
        f << std::setprecision(precision) << A[i] << ",";
      else
        f << "0"
          << ",";
      f << std::setprecision(precision) << B[i] << ",";
      if (i != N - 1)
        f << std::setprecision(precision) << C[i] << std::endl;
      else
        f << 0 << std::endl;
    }
    f.close();
  }

  // Solves the system Ax = d in place where d eventually stores the solution
  void solveTDMA(std::vector<T> & d)
  {
    // Forward sweep
    T tmp = 0;
    for (unsigned int i = 1; i < N; ++i)
    {
      tmp = A[i] / B[i - 1];
      B[i] -= tmp * C[i - 1];
      d[i] -= tmp * d[i - 1];
    }

    // Backward sweep
    d[N - 1] /= B[N - 1];
    for (unsigned int i = N - 2; i != std::numeric_limits<unsigned int>::max(); --i)
    {
      d[i] -= C[i] * d[i + 1];
      d[i] /= B[i];
    }
  }

protected:
  // Matrix size (N x N)
  unsigned int N;

  // Left/main/right diagonal storage
  std::vector<T> A, B, C;
};

#endif /* TRIDIAGONAL_H */
