// Evaluate a polynomial.
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "utils.h"

double polyeval(Eigen::VectorXd coeffs, double x) {
    double result = 0.0;
    for (int i = 0; i < coeffs.size(); i++) {
      result += coeffs[i] * pow(x, i);
    }
    return result;
  }
  
  // Fit a polynomial.
  // Adapted from
  // https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
  Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                          int order) {
    assert(xvals.size() == yvals.size());
    assert(order >= 1 && order <= xvals.size() - 1);
    Eigen::MatrixXd A(xvals.size(), order + 1);
  
    for (int i = 0; i < xvals.size(); i++) {
      A(i, 0) = 1.0;
    }
  
    for (int j = 0; j < xvals.size(); j++) {
      for (int i = 0; i < order; i++) {
        A(j, i + 1) = A(j, i) * xvals(j);
      }
    }
  
    auto Q = A.householderQr();
    auto result = Q.solve(yvals);
    return result;
  }
  
  double find_slope(Eigen::VectorXd curve, double x){
    Eigen::VectorXd coeff_derivative(curve.size()-1);
    for(int i =0; i<coeff_derivative.size(); i++){
      coeff_derivative(i) = (i+1)*curve(i+1);
    }
  
    double slope = polyeval(coeff_derivative,x);
    return slope;
  }