#include "Base.h"
#include "Conduction2D.h"

int main() {
  Conduction2D problem(50, 50, 1.0);

  problem.run();
  problem.save("test.csv");

  return 0;
}
