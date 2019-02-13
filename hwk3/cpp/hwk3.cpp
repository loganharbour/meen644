#include "Conduction2D.h"

int main() {
  // Problem 1
  for (double w : {1.0, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40})
  {
    std::cout << "Running with w = " << w << ": ";
    Conduction2D problem(15, 15, w);
    problem.run();
  }
  std::cout << std::endl;

  // Problem 2
  for (unsigned int N : {15, 21, 25, 31, 41}){
    std::cout << "Running with " << N * N << " CVs: ";
    auto center = (N - 1) / 2;
    Conduction2D problem(N, N, 1.2);
    problem.run();
    std::cout << "  Center solution = " << std::setprecision(6)
              << problem.getT()(center, center) << " C" << std::endl;

    // Problem 3
    if (N == 41)
      problem.getT().save("solution.csv");
  }

  return 0;
}
