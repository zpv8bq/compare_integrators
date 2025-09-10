#include <stdlib.h>
#include <cmath>
#include <stdio.h>
#include "integrators.hpp"
#include "gaussrules.hpp"
#include <iostream>
#include <iomanip>
using std::cout;
using std::endl;

// Note: this code is not safe wrt the number of intervals
// Make sure the number of intervals used is valid for each integrator


// Integration using trapezoid rule 
// npts : number of points used in calculation (npts>=2)
double trapez (double (*f)(double x), unsigned npts, double min, double max) {
  double sum=0.;		 

  // complete your code here
  
  return (sum);
}      

// Integration using Simpson's rule
// npts : number of points used in calculation (npts odd, and >=3)
double simpson (double (*f)(double x), unsigned npts, double min, double max){  
  double sum=0.;

  // complete your code here
  
  return (sum);
}  

// Integration using Gauss's rule, code is based on the Landau text
// This is not a very good implementation
double gaussint (double (*f)(double x), unsigned npts, double min, double max){
  double result = 0.;
  const unsigned MAXPOINTS=10000;
  double w[MAXPOINTS];    // for points and weights
  double xi[MAXPOINTS];
  
  if (npts>MAXPOINTS) npts=MAXPOINTS;
  gauss (npts, min, max, xi, w);      // returns Legendre polynomials
  // points and weights
  double c1 = (max - min) / 2;
  double c2 = (max + min) / 2;
  for (unsigned n=0; n<npts; n++){
    result += w[n] * f(c1 * xi[n] + c2);   // calculating the integral
  } 
  return result/2;                  
}

// calculate the weights and intervals for the n-point rule
void gauss(unsigned npts, double a, double b, double x[], double w[]){    
  // npts     number of points
  // x, w     output grid interval points and weights.			      

  double  t, t1, pp=0, p1, p2, p3;
  double  eps = 3.e-10;			// limit for accuracy

  // calculating roots of Legendre polynomials and wgt cofficients
  unsigned m = (npts+1)/2;
  for(unsigned i=1; i<=m; i++){  
    t  = cos(M_PI*(i-0.25)/(npts+0.5));
    t1 = 1;
    while(fabs(t-t1)>=eps){ 
      p1 = 1.0;
      p2 = 0.0;
      for(unsigned j=1; j<=npts; j++) {
	p3 = p2;
	p2 = p1;
	p1 = ((2*j-1)*t*p2-(j-1)*p3)/j;
      }
      pp = npts*(t*p1-p2)/(t*t-1);
      t1 = t;
      t  = t1 - p1/pp;
    }   
    x[i-1] = -t;
    x[npts-i] = t;
    w[i-1]    = 2.0/((1-t*t)*pp*pp);
    w[npts-i] = w[i-1];
  } 
}

// modified from rosettacode.org
// this is a MUCH better implemtation of the Gaussian quadrature method
GaussInt::GaussInt(int npoints){
  Init(npoints);
}

#if __GNUC_PREREQ(4,6)
__float128 GaussInt::lege_eval(int n, __float128 x){
  __float128 s = lcoef[n][n];
  for (int i = n; i>0; i--)
    s = s * x + lcoef[n][i - 1];
  return s;
}

__float128 GaussInt::lege_diff(int n, __float128 x){
  return n * (x * lege_eval(n, x) - lege_eval(n - 1, x)) / (x * x - 1);
}
#endif


void GaussInt::Init(int npoints){

  // calculates abscissas and weights to double precision accuracy
  // for n-point quadrature rule
  // Note: In practice you would generally not bother with this step,
  // instead tables of abscissas and weights would be calculated once for
  // a selection of n-point rules and saved for later usage

#if __GNUC_PREREQ(4,6)  
  lroots.assign(npoints,0);
  weight.assign(npoints,0);
  lcoef.assign(npoints+1,vector<__float128>(npoints+1,0));
  
  lcoef[0][0] = lcoef[1][1] = 1;
  for (int n = 2; n <= npoints; n++) {
    lcoef[n][0] = -(n - 1) * lcoef[n - 2][0] / n;
    for (int i = 1; i <= n; i++)
      lcoef[n][i] = ((2 * n - 1) * lcoef[n - 1][i - 1]
		     - (n - 1) * lcoef[n - 2][i] ) / n;
  }

  __float128 x, x1;
  for (int i = 1; i <= npoints; i++) {
    x = cos(M_PI * (i - .25) / (npoints + .5));
    do {
      x1 = x;
      x -= lege_eval(npoints, x) / lege_diff(npoints, x);
    } while (fabs((long double)(x-x1))>1e-16);  // keep double precision 
    lroots[i - 1] = x;
    
    x1 = lege_diff(npoints, x);
    weight[i - 1] = 2 / ((1 - x * x) * x1 * x1);
  }
  // cout << "==== " << npoints << " ====" << endl;
  // cout.precision(20);
  // for (int i = 0; i < npoints; i++)
  //   cout << (double) weight[i] << " " << (double) lroots[i] << endl;
#else
  // if we have an older compiler use the precomputed values
  if (npoints>MAXPOINTRULE){
    npoints=MAXPOINTRULE;
    cout << "Truncating to " << MAXPOINTRULE << "-point rule" << endl;
    cout << "Use gcc 4.6 or higher to enable calculation of higher point rules" << endl;
  }
  double *wa = &GAUSSWA[npoints*npoints-(npoints+2)];
  weight.clear();
  lroots.clear();
  for (int i=0; i<npoints*2; i+=2){
    weight.push_back(wa[i]);
    lroots.push_back(wa[i+1]);
  }
#endif  
}

// notice the efficiency of the integration
// once the constants for the n-point rule are calculated
// the integral is a very simple sum
double GaussInt::Integ(double (*f)(double x), double a, double b){

  // ### complete code here ###
  
  return 0;  // integral of function f
}

void GaussInt::PrintWA() const{
  cout << "Weights, abscissas for " << weight.size() << " point rule" << endl;
  for (unsigned i = 0; i < weight.size(); i++){
    cout  << std::setprecision(34) << (long double) weight[i] << "   " << (long double) lroots[i] << endl;
  }
}
