#include "Conduction2D.h"
#include <boost/format.hpp>
#include <map>
#include <sstream>

int main() {
  // Part a: change relaxation factor
  std::cout << "Part a" << std::endl;
  for (double w : {1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4}) {
    std::cout << "  Running with w = " << w << ": ";
    Conduction2D problem(15, 15, w);
    problem.solve();

    // Save residuals if converged
    if (problem.converged()) {
      std::stringstream filename;
      filename << boost::format("results/a/w%1%.csv") % w;
      saveCSV(problem.getResiduals(), filename.str());
    }
  }
  std::cout << std::endl;

  // Part b: refine mesh and retreive center value
  std::cout << "Part b" << std::endl;
  std::map<double, std::vector<double>> temps;
  for (unsigned int N : {15, 21, 25, 31, 41}) {
    std::cout << "  Running with " << N * N << " CVs: " << std::endl;
    for (double w : {1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3}) {
      std::cout << "    Running with w = " << w << ": ";

      Conduction2D problem(N, N, w);
      problem.solve();

      // Store iteration count and center solution if it converged
      if (problem.converged()) {
        std::stringstream filename;
        filename << boost::format("results/b-iterations/N%1%_w%2%.csv") % N % w;
        saveCSV(problem.getResiduals(), filename.str());
        temps[w].push_back(problem.getT((N - 1) / 2, (N - 1) / 2));
      }

      // Part c: save solution for 41x41
      if (N == 41 && w == 1.0)
        problem.saveT("results/c/solution.csv");
    }
  }
  // Save temperatures for each alpha
  for (double w : {1.0, 1.15, 1.3})
    saveCSV(temps[w], "results/b-temps/" + std::to_string(w) + ".csv");

  return 0;
}
