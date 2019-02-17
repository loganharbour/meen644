#ifndef MATRIX
#define MATRIX

#define NDEBUG
#include <cassert>
#include <vector>

/**
 * Class that holds a Nx x Ny matrix with common matrix operations.
 */
class Matrix {
public:
  Matrix(unsigned int Nx, unsigned int Ny, double v = 0)
      : Nx(Nx), Ny(Ny), M(Nx, std::vector<double>(Ny, v)) {}

  // Const operator for getting the (i, j) element
  const double &operator()(unsigned int i, unsigned int j) const {
    assert(i < Nx && j < Ny);
    return M[i][j];
  }
  // Operator for getting the (i, j) element
  double &operator()(unsigned int i, unsigned int j) {
    assert(i < Nx && j < Ny);
    return M[i][j];
  }
  // Operator for setting the entire matrix to a value
  void operator=(double v) {
    for (unsigned int j = 0; j < Ny; ++j)
      setRow(j, v);
  }

  // Saves the matrix in csv format
  void save(const std::string filename, unsigned int precision = 12) const {
    std::ofstream f;
    f.open(filename);
    for (unsigned int j = 0; j < Ny; ++j) {
      for (unsigned int i = 0; i < Nx; ++i) {
        if (i > 0)
          f << ",";
        f << std::setprecision(precision) << M[i][j];
      }
      f << std::endl;
    }
    f.close();
  }

  // Set the j-th row to v
  void setRow(unsigned int j, double v) {
    assert(j < Ny);
    for (unsigned int i = 0; i < Nx; ++i)
      M[i][j] = v;
  }
  // Set the i-th column to v
  void setColumn(unsigned int i, double v) {
    assert(i < Nx);
    for (unsigned int j = 0; j < Ny; ++j)
      M[i][j] = v;
  }

  // Set the j-th row to vs
  void setRow(unsigned int j, std::vector<double> &vs) {
    assert(j < Ny && vs.size() == Nx);
    for (unsigned int i = 0; i < Nx; ++i)
      M[i][j] = vs[i];
  }
  // Set the i-th column to vs
  void setColumn(unsigned int i, std::vector<double> &vs) {
    assert(i < Nx && vs.size() == Ny);
    for (unsigned int j = 0; j < Ny; ++j)
      M[i][j] = vs[j];
  }

private:
  // The size of this matrix
  const unsigned int Nx, Ny;

  // Matrix storage
  std::vector<std::vector<double> > M;
};

#endif /* MATRIX_H */
