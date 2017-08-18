#include "opt.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"
#include <iostream>

using CppAD::AD;

// Set the timestep length
double dt = 0.02;

// This value assumes the model presented in the classroom is used.
//
// It was obtained by measuring the radius formed by running the vehicle in the
// simulator around in a circle with a constant steering angle and velocity on a
// flat terrain.
//
// Lf was tuned until the the radius formed by the simulating the model
// presented in the classroom matched the previous radius.
//
// This is the length from front to CoG that has a similar radius.
const double Lf = 2.67;

// Both the reference cross track and orientation errors are 0.
// The reference velocity is set to 40 mph.
double ref_v = 49.5 * 0.44704;
double ref_d = 2.0;

// The solver takes all the state variables and actuator
// variables in a singular vector. Thus, we should to establish
// when one variable starts and another ends to make our lifes easier.
#define s_start 0
#define d_start s_start + N
#define vs_start d_start + N
#define vd_start vs_start + N
#define as_start vd_start + N
#define ad_start as_start + N
#define js_start ad_start + N
#define jd_start js_start + N

class FG_eval {
 public:
  // Fitted polynomial paramicients
  Eigen::VectorXd params;
  tk::spline waypoint_spline;
  FG_eval(Eigen::VectorXd params){
    this->params = params;
  }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // TODO: implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.

    int N = params[0];
    double ref_d = params[1];
    // The cost is stored is the first element of `fg`.
    // Any additions to the cost should be added to `fg[0]`.
    fg[0] = 0;

    // The part of the cost based on the reference state.
    for (int t = 0; t < N; t++) {
      // Add cost to accelerations
      fg[0] += CppAD::pow(vars[as_start + t], 2);
      fg[0] += CppAD::pow(vars[ad_start + t], 2);
      // Add cost to jerks
      fg[0] += CppAD::pow(vars[js_start + t], 2);
      fg[0] += CppAD::pow(vars[jd_start + t], 2);
      // Add cost to difference between reference velocity in s_direction.
      fg[0] += 10 * CppAD::pow(
        CppAD::sqrt(CppAD::pow(vars[vs_start + t], 2)+CppAD::pow(vars[vd_start + t], 2))
         - ref_v,2);
      // Add cost to difference between reference lane positon in d_direction.
      fg[0] += 10 * CppAD::pow(vars[d_start + t] - ref_d, 2);
    }

    //
    // Setup Constraints
    //
    // NOTE: In this section you'll setup the model constraints.

    fg[1 + s_start] = vars[s_start];
    fg[1 + d_start] = vars[d_start];
    fg[1 + vs_start] = vars[vs_start];
    fg[1 + vd_start] = vars[vd_start];
    fg[1 + as_start] = vars[as_start];
    fg[1 + ad_start] = vars[ad_start];

    // Integral Constraints
    for (int t = 1; t < N; t++) {
      // The state at time t+1 .
      AD<double> s1 = vars[s_start + t];
      AD<double> d1 = vars[d_start + t];
      AD<double> vs1 = vars[vs_start + t];
      AD<double> vd1 = vars[vd_start + t];
      AD<double> as1 = vars[as_start + t];
      AD<double> ad1 = vars[ad_start + t];
      AD<double> js1 = vars[js_start + t];
      AD<double> jd1 = vars[jd_start + t];

      // The state at time t.
      AD<double> s0 = vars[s_start + t - 1];
      AD<double> d0 = vars[d_start + t - 1];
      AD<double> vs0 = vars[vs_start + t - 1];
      AD<double> vd0 = vars[vd_start + t - 1];
      AD<double> as0 = vars[as_start + t - 1];
      AD<double> ad0 = vars[ad_start + t - 1];
      AD<double> js0 = vars[js_start + t - 1];
      AD<double> jd0 = vars[jd_start + t - 1];

      // Integral constraints
      fg[1 + s_start + t] = s1 - (s0 + vs0 * dt);
      fg[1 + d_start + t] = d1 - (d0 + vd0 * dt);
      fg[1 + vs_start + t] = vs1 - (vs0 + as0 * dt);
      fg[1 + vd_start + t] = vd1 - (vd0 + ad0 * dt);
      fg[1 + as_start + t] = as1 - (as0 + js0 * dt);
      fg[1 + ad_start + t] = ad1 - (ad0 + jd0 * dt);

    }
  }
};

//
// MPC class definition implementation.
//
OPT::OPT() {is_initial = true;}
OPT::~OPT() {}

vector<double> OPT::Solve(Eigen::VectorXd state, Eigen::VectorXd params) {
  bool ok = true;
  size_t i;
  typedef CPPAD_TESTVECTOR(double) Dvector;
  this->is_initial = false;

  double s = state[0];
  double d = state[1];
  double vs = state[2];
  double vd = state[3];
  double as = state[4];
  double ad = state[5];

  int N = params[0];

  // Number of variable
  size_t n_vars = N * 8;
  // Number of constraints for g(x)
  size_t n_constraints = N * 6;


  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (int i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }

  // Set the initial variable values
  vars[s_start] = s;
  vars[d_start] = d;
  vars[vs_start] = vs;
  vars[vd_start] = vd;
  vars[as_start] = as;
  vars[ad_start] = ad;


  // Set lower and upper boundary for variables
  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);

  for (int t = 0; t < N; t++){
    vars_lowerbound[s_start + t] = -1.0e19;
    vars_lowerbound[d_start + t] = -1.0e19;
    vars_lowerbound[vs_start + t] = -1.0e19;
    vars_lowerbound[vd_start + t] = -1.0e19;
    vars_lowerbound[as_start + t] = -5.0;
    vars_lowerbound[ad_start + t] = -5.0;
    vars_lowerbound[js_start + t] = -5.0;
    vars_lowerbound[jd_start + t] = -5.0;

    vars_upperbound[s_start + t] = 1.0e19;
    vars_upperbound[d_start + t] = 1.0e19;
    vars_upperbound[vs_start + t] = 1.0e19;
    vars_upperbound[vd_start + t] = 1.0e19;
    vars_upperbound[as_start + t] = 5.0;
    vars_upperbound[ad_start + t] = 5.0;
    vars_upperbound[js_start + t] = 5.0;
    vars_upperbound[jd_start + t] = 5.0;
  }


  // Set lower and upper boundary for the constraints
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (int i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }
  constraints_lowerbound[s_start] = s;
  constraints_lowerbound[d_start] = d;
  constraints_lowerbound[vs_start] = vs;
  constraints_lowerbound[vd_start] = vd;
  constraints_lowerbound[as_start] = as;
  constraints_lowerbound[ad_start] = ad;

  constraints_upperbound[s_start] = s;
  constraints_upperbound[d_start] = d;
  constraints_upperbound[vs_start] = vs;
  constraints_upperbound[vd_start] = vd;
  constraints_upperbound[as_start] = as;
  constraints_upperbound[ad_start] = ad;

  // object that computes objective and constraints
  FG_eval fg_eval(params);

  // options for IPOPT solver
  std::string options;
  // Uncomment this if you'd like more print information
  options += "Integer print_level  0\n";
  // NOTE: Setting sparse to true allows the solver to take advantage
  // of sparse routines, this makes the computation MUCH FASTER. If you
  // can uncomment 1 of these and see if it makes a difference or not but
  // if you uncomment both the computation time should go up in orders of
  // magnitude.
  options += "Sparse  true        forward\n";
  options += "Sparse  true        reverse\n";
  // NOTE: Currently the solver has a maximum time limit of 0.5 seconds.
  // Change this as you see fit.
  options += "Numeric max_cpu_time          0.5\n";

  // place to return solution
  CppAD::ipopt::solve_result<Dvector> solution;

  // solve the problem
  cout << "Start solving" << endl;
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);
  cout << "Solved" << endl << endl;

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;

  // TODO: Return the first actuator values. The variables can be accessed with
  // `solution.x[i]`.
  //
  // {...} is shorthand for creating a vector, so auto x1 = {1.0,2.0}
  // creates a 2 element double vector.
  this->solution.clear();
  this->solution.shrink_to_fit();
  for(int i=0; i<N*8; ++i){
    this->solution.push_back(solution.x[i]); // IPOPT solver result
  }

  return this->solution;
}
