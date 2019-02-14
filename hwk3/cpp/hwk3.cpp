#include "Conduction2D.h"

int main() {
  // Part a: change relaxation factor
  std::cout << "Part a" << std::endl;
  for (double w : {1.0, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40}) {
    std::cout << " Running with w = " << w << ": ";
    Conduction2D problem(15, 15, w);
    problem.solve();
  }
  std::cout << std::endl;

  // Part b: refine mesh and retreive center value
  std::cout << "Part b" << std::endl;
  for (unsigned int N : {15, 21, 25, 31, 41}) {
    std::cout << " Running with " << N * N << " CVs: ";

    Conduction2D problem(N, N, 1.2);
    problem.solve();
    std::cout << "  Center solution = " << std::setprecision(10)
              << problem.getT()((N - 1) / 2, (N - 1) / 2) << " C" << std::endl;

    // Part c: save finest solution
    if (N == 41)
      problem.save("solution.csv");
  }

  return 0;
}
