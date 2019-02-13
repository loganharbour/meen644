#ifndef MATRIX
#define MATRIX

#define NDEBUG

#include <vector>

class Matrix {
public:
  Matrix(unsigned int Nx, unsigned int Ny, double val = 0)
      : Nx(Nx), Ny(Ny), M(Nx, std::vector<double>(Ny, val)) {}

  const double &operator()(unsigned int i, unsigned int j) const {
    assert(i < Nx && j < Ny);
    return M[i][j];
  }
  double &operator()(unsigned int i, unsigned int j) {
    assert(i < Nx && j < Ny);
    return M[i][j];
  }
  Matrix &operator=(const double rhs) {
    setAll(rhs);
    return *this;
  }

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

  void setRow(unsigned int j, double val) {
    assert(j < Ny);
    for (unsigned int i = 0; i < Nx; ++i)
      M[i][j] = val;
  }
  void setColumn(unsigned int i, double val) {
    assert(i < Nx);
    for (unsigned int j = 0; j < Ny; ++j)
      M[i][j] = val;
  }
  void setAll(double val) {
    for (unsigned int j = 0; j < Ny; ++j)
      setRow(j, val);
  }

  void setRow(unsigned int j, std::vector<double> &vals) {
    assert(j < Ny && vals.size() == Nx);
    for (unsigned int i = 0; i < Nx; ++i)
      M[i][j] = vals[i];
  }
  void setColumn(unsigned int i, std::vector<double> &vals) {
    assert(i < Nx && vals.size() == Ny);
    for (unsigned int j = 0; j < Ny; ++j)
      M[i][j] = vals[j];
  }

private:
  const unsigned int Nx, Ny;
  std::vector<std::vector<double> > M;
};

#endif /* MATRIX_H */
