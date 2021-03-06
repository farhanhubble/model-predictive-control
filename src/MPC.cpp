#include "MPC.h"
#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"
#include "utils.h"

using CppAD::AD;

AD<double> polyeval_AD(Eigen::VectorXd coeffs, AD<double> x){
  AD<double> result = 0;
  for (int i = 0; i < coeffs.size(); i++) {
    result += coeffs[i] * CppAD::pow(x, i);
  }
  return result;
}


AD<double>  find_slope_AD(Eigen::VectorXd curve, AD<double>  x){
  Eigen::VectorXd coeff_derivative(curve.size()-1);
  for(int i =0; i<coeff_derivative.size(); i++){
    coeff_derivative(i) = (i+1)*curve(i+1);
  }

  AD<double>  slope = polyeval_AD(coeff_derivative,x);
  return slope;
}

// TODO: Set the timestep length and duration
size_t N = 15;
double dt = 0.1;

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


// Reference velocity.
const double v_ref = 40.0;

size_t x_start = 0;
size_t y_start = x_start + N;
size_t psi_start = y_start + N;
size_t v_start = psi_start + N;
size_t cte_start = v_start + N;
size_t epsi_start = cte_start + N;
size_t delta_start = epsi_start + N;
size_t a_start = delta_start + N - 1;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;
  FG_eval(Eigen::VectorXd coeffs) { this->coeffs = coeffs; }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector& fg, const ADvector& vars) {
    // Initialize cost.
    fg[0] = 0;
    
    // Add deviation from reference state to cost.
    for (size_t t = 0; t < N; t++) {
      fg[0] += 20*CppAD::pow(vars[cte_start + t], 2);
      fg[0] += 20*CppAD::pow(vars[epsi_start + t], 2);
      fg[0] += 10*CppAD::pow(vars[v_start + t] - v_ref, 2);
    }

    // Penalize the use of actuators.
    for (size_t t = 0; t < N - 1; t++) {
      fg[0] += 1000*CppAD::pow(vars[delta_start + t], 2);
      fg[0] += 1000*CppAD::pow(vars[a_start + t], 2);
    }

    // Penalize sudden changes in sequential actuations.
    for (size_t t = 0; t < N - 2; t++) {
      fg[0] += 500*CppAD::pow(vars[delta_start + t + 1] - vars[delta_start + t], 2);
      fg[0] += 500*CppAD::pow(vars[a_start + t + 1] - vars[a_start + t], 2);
    }

    // Set initial constraints.
    fg[1 + x_start] = vars[x_start];
    fg[1 + y_start] = vars[y_start];
    fg[1 + psi_start] = vars[psi_start];
    fg[1 + v_start] = vars[v_start];
    fg[1 + cte_start] = vars[cte_start];
    fg[1 + epsi_start] = vars[epsi_start];

    // The rest of the constraints.
    // The constraints are calculated as the difference between the predicted 
    // and 
    for (size_t t = 1; t < N; t++) {

      // The state at time t+1 .
      AD<double> x1 = vars[x_start + t];
      AD<double> y1 = vars[y_start + t];
      AD<double> psi1 = vars[psi_start + t];
      AD<double> v1 = vars[v_start + t];
      AD<double> cte1 = vars[cte_start + t];
      AD<double> epsi1 = vars[epsi_start + t];

      // The state at time t.
      AD<double> x0 = vars[x_start + t - 1];
      AD<double> y0 = vars[y_start + t - 1];
      AD<double> psi0 = vars[psi_start + t - 1];
      AD<double> v0 = vars[v_start + t - 1];
      AD<double> cte0 = vars[cte_start + t - 1];
      AD<double> epsi0 = vars[epsi_start + t - 1];

      // Only consider the actuation at time t.
      AD<double> delta0 = vars[delta_start + t - 1];
      AD<double> a0 = vars[a_start + t - 1];
      AD<double> f0 = polyeval_AD(coeffs,x0);
      AD<double> psides0 = CppAD::atan(find_slope_AD(coeffs,x0));


      fg[1 + x_start + t]    = x1 - (x0 + v0 * CppAD::cos(psi0) * dt);
      fg[1 + y_start + t]    = y1 - (y0 + v0 * CppAD::sin(psi0) * dt);
      fg[1 + psi_start + t]  = psi1 - (psi0 + v0 * delta0 / Lf * dt);
      fg[1 + v_start + t]    = v1 - (v0 + a0 * dt);
      fg[1 + cte_start + t]  = cte1 - ((f0 - y0) + (v0 * CppAD::sin(epsi0) * dt));
      fg[1 + epsi_start + t] = epsi1 - ((psi0 - psides0) + v0 * delta0 / Lf * dt);
    }
  }
};

//
// MPC class definition implementation.
//
MPC::MPC() {}
MPC::~MPC() {}

vector<double> MPC::Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs) {
  bool ok = true;
  typedef CPPAD_TESTVECTOR(double) Dvector;

  // Set the total number of model variables. 
  // Total model variables = #state variables * #timesteps + #actuation variables * (#timesteps -1)  
  const size_t n_actuators = 2;
  size_t n_vars = state.size() * N + n_actuators * (N-1);

  // Set total number of constraints. 
  // There is one constraint per state variable per time step.
  size_t n_constraints = state.size() * N ;

  // Extract initial state.
  double x    = state[0];
  double y    = state[1];
  double psi  = state[2];
  double v    = state[3];
  double cte  = state[4];
  double epsi = state[5];

  // Initial value of the independent variables.
  // SHOULD BE 0 besides initial state.
  Dvector vars(n_vars);
  for (size_t i = 0; i < n_vars; i++) {
    vars[i] = 0;
  }
  vars[x_start] = x;
  vars[y_start] = y;
  vars[psi_start] = psi;
  vars[v_start] = v;
  vars[cte_start] = cte;
  vars[epsi_start] = epsi;

  //Define upper and lower limits for all variables.

  // Set all non-actuator upper and lower limits
  // to the max negative and positive values.
  Dvector vars_lowerbound(n_vars);
  Dvector vars_upperbound(n_vars);
  for (size_t i = 0; i < delta_start; i++) {
    vars_lowerbound[i] = -1.0e19;
    vars_upperbound[i] = 1.0e19;
  }

  // The upper and lower limits of steering angle are set to -25 and 25
  // degrees (converted to radians).
  for (size_t i = delta_start; i < a_start; i++) {
    vars_lowerbound[i] = -0.436332;
    vars_upperbound[i] = 0.436332;
  }

  // Throttle upper and lower limits.
  for (size_t i = a_start; i < n_vars; i++) {
    vars_lowerbound[i] = -1.0;
    vars_upperbound[i] = 1.0;
  }

  // Lower and upper limits for the constraints
  // Should be 0 besides initial state.
  Dvector constraints_lowerbound(n_constraints);
  Dvector constraints_upperbound(n_constraints);
  for (size_t i = 0; i < n_constraints; i++) {
    constraints_lowerbound[i] = 0;
    constraints_upperbound[i] = 0;
  }

  constraints_lowerbound[x_start] = x;
  constraints_lowerbound[y_start] = y;
  constraints_lowerbound[psi_start] = psi;
  constraints_lowerbound[v_start] = v;
  constraints_lowerbound[cte_start] = cte;
  constraints_lowerbound[epsi_start] = epsi;

  constraints_upperbound[x_start] = x;
  constraints_upperbound[y_start] = y;
  constraints_upperbound[psi_start] = psi;
  constraints_upperbound[v_start] = v;
  constraints_upperbound[cte_start] = cte;
  constraints_upperbound[epsi_start] = epsi;

  // object that computes objective and constraints
  FG_eval fg_eval(coeffs);


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
  CppAD::ipopt::solve<Dvector, FG_eval>(
      options, vars, vars_lowerbound, vars_upperbound, constraints_lowerbound,
      constraints_upperbound, fg_eval, solution);

  // Check some of the solution values
  ok &= solution.status == CppAD::ipopt::solve_result<Dvector>::success;

  // Cost
  auto cost = solution.obj_value;
  std::cout << "Cost " << cost << std::endl;
  vector<double> ret;
  
    ret.push_back(solution.x[delta_start]);
    ret.push_back(solution.x[a_start]);
  
    for (size_t i = 0; i < N-1; i++) {
      ret.push_back(solution.x[1 + x_start + i ]);
      ret.push_back(solution.x[1 + y_start + i ]);
    }
  
  return ret;
}
