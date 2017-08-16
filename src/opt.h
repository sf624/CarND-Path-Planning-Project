#ifndef OPT_H
#define OPT_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

class OPT {
 public:
  OPT();
  bool is_initial;
  vector<double> solution;

  virtual ~OPT();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd params);
};

#endif /* OPT_H */
