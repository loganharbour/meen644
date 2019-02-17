#include "Conduction2D.h"
#include <map>

template <typename T> void save(const std::vector<T> &v, std::string filename) {
  std::ofstream f;
  f.open(filename);
  for (unsigned int i = 0; i < v.size(); ++i)
    f << std::scientific << v[i] << std::endl;
  f.close();
}

int main() {
  // Part a: change relaxation factor
  std::cout << "Part a" << std::endl;
  for (double w : {1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4}) {
    std::cout << "  Running with w = " << w << ": ";
    Conduction2D problem(15, 15, w);
    problem.solve();
    save(problem.getResiduals(), "results/a/" + std::to_string(w) + ".csv");
  }
  std::cout << std::endl;

  // Part b: refine mesh and retreive center value
  std::cout << "Part b" << std::endl;
  std::map<double, std::vector<double>> temps;
  for (unsigned int N : {15, 21, 25, 31, 41}) {
    std::cout << "  Running with " << N * N << " CVs: " << std::endl;

    std::vector<unsigned int> iterations;

    for (double w : {1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35}) {
      std::cout << "    Running with w = " << w << ": ";

      Conduction2D problem(N, N, w);
      problem.solve();

      // Store iteration count and center solution if it converged
      if (problem.getNumIterations() != 1000) {
        iterations.push_back(problem.getNumIterations());
        temps[w].push_back(problem.getT((N - 1) / 2, (N - 1) / 2));
      }

      // Part c: save solution for 41x41
      if (N == 41 && w == 1.0)
        problem.saveT("results/c/solution.csv");
    }

    // Save iteration count and center solution
    save(iterations, "results/b-iterations/" + std::to_string(N) + ".csv");
  }
  for (double w : {1.0, 1.15, 1.3})
    save(temps[w], "results/b-temps/" + std::to_string(w) + ".csv");
  return 0;
}
