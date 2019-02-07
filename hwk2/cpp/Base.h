#ifndef BASE_H
#define BASE_H

#define NDEBUG

#include <cassert>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>

/**
 * Saves a vector to a csv.
 */
void saveVectorCsv(const std::vector<double> x, const std::string filename) {
  std::ofstream f;
  f.open(filename);
  for (unsigned int i = 0; i < x.size(); ++i)
    f << std::setprecision(12) << x[i] << std::endl;
  f.close();
}

#endif /* BASE_H */
