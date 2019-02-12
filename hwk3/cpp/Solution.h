#ifndef SOLUTION_H
#define SOLUTION_H

#include "Base.h"

class Solution {
public:
  Solution(unsigned int Nx, unsigned int Ny)
      : _sol(Nx * Ny), _Nx(Nx), _Ny(Ny) {}

  const double get(unsigned int i, unsigned int j) {
    assert(i < _Nx && j < _Ny);
    return _sol[j * _Nx + i];
  }
  const void getRow(unsigned int row, std::vector<double> &vals) {
    assert(row < _Ny && vals.size() == _Nx);
    std::copy_n(_sol.begin() + row * _Nx, _Nx, vals.begin());
  }
  const void getColumn(unsigned int col, std::vector<double> &vals) {
    assert(col < _Nx && vals.size() == _Ny);
    for (unsigned int i = 0; i < _Ny; ++i)
      vals[i] = _sol[col + _Nx * i];
  }

  // Setters for rows and columns
  void set(unsigned int i, unsigned int j, double val) {
    assert(i < _Nx && j < _Ny);
    _sol[j * _Nx + i] = val;
  }
  void setRow(unsigned int row, std::vector<double> &vals) {
    assert(row < _Ny && vals.size() == _Nx);
    std::copy_n(vals.begin(), _Nx, _sol.begin() + row * _Nx);
  }
  void setColumn(unsigned int col, std::vector<double> &vals) {
    assert(col < _Nx && vals.size() == _Ny);
    for (unsigned int i = 0; i < _Ny; ++i)
      _sol[col + _Nx * i] = vals[i];
  }

  const std::vector<double> &get() { return _sol; }
  const void save(std::string filename) { saveVectorCsv(_sol, filename); }

private:
  unsigned int _Nx, _Ny;
  std::vector<double> _sol;
};

#endif /* SOLUTION_H */
