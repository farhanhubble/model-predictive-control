#ifndef __UTILS_H__
#define __UTILS_H__
// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x);
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,int order);
double find_slope(Eigen::VectorXd curve, double x);

#endif  // __UTILS_H__