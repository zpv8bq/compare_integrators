#ifndef INTEGRATORS_HPP
#define INTEGRATORS_HPP

#include <vector>
#include <features.h>
using std::vector;

// Integration using trapezoid rule
double trapez (double (*f)(double x), unsigned npts, double min, double max);

// Integration using Simpson's rule 
double simpson (double (*f)(double x), unsigned npts, double min, double max);

// Integration using Gauss's rule 
double gaussint (double (*f)(double x), unsigned npts, double min, double max);

// calculation of abscissas and weights
// not a very good implementation from the text!
void gauss(unsigned npts, double a, double b, double x[], double w[]);


// this version does a MUCH better job of determining 
// abscissas and weights.  Of course these can (and should) be tabulated 1 time
// and stored as constants in a header file!
class GaussInt {
public:
  GaussInt(int npoints=5);
  void Init(int npoints);
  double Integ(double (*f)(double x), double a, double b);
  void PrintWA() const;
private:
#if __GNUC_PREREQ(4,6) // good calculation requires support for quad precision
  __float128 lege_eval(int n, __float128 x);
  __float128 lege_diff(int n, __float128 x);
  vector<__float128> lroots;
  vector<__float128> weight;
  vector<vector<__float128> > lcoef;
#else
  vector<double> lroots;
  vector<double> weight;
#endif
};
#endif


