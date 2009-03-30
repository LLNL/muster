/** \file ltqnorm.C
 * 
 * Based on the code by Peter J. Acklam, and the adaptation by V. Natarajan,
 * available from http://home.online.no/~pjacklam/notes/invnorm/.
 * The code in this file is free to distribute or use for any purpose.
 *
 * Adapted for C++ by Todd Gamblin, 1/28/2005.  tgamblin@cs.unc.edu
 *
 * Original comments follow --------------------------------------------
 *
 * Z = LTQNORM(P) returns the lower tail quantile for the standard normal
 * distribution function.  I.e., it returns the Z satisfying Pr{X < Z} = P,
 * where X has a standard normal distribution.
 *
 * LTQNORM(P) is the same as SQRT(2) * ERFINV(2*P-1), but the former returns a
 * more accurate value when P is close to zero.
 *
 * The algorithm uses a minimax approximation by rational functions and the
 * result has a relative error less than 1.15e-9.  A last refinement by
 * Halley's rational method is applied to achieve full machine precision.
 *
 * Author:      Peter J. Acklam
 * Time-stamp:  2003-04-23 08:26:51 +0200
 * E-mail:      pjacklam@online.no
 * URL:         http://home.online.no/~pjacklam
 */
#include <cmath>
#include "ltqnorm.h"


//approximation constants for different parts of the function
static const double  A1 = -3.969683028665376e+01;
static const double  A2 =  2.209460984245205e+02;
static const double  A3 = -2.759285104469687e+02;
static const double  A4 =  1.383577518672690e+02;
static const double  A5 = -3.066479806614716e+01;
static const double  A6 =  2.506628277459239e+00;

static const double  B1 = -5.447609879822406e+01;
static const double  B2 =  1.615858368580409e+02;
static const double  B3 = -1.556989798598866e+02;
static const double  B4 =  6.680131188771972e+01;
static const double  B5 = -1.328068155288572e+01;

static const double  C1 = -7.784894002430293e-03;
static const double  C2 = -3.223964580411365e-01;
static const double  C3 = -2.400758277161838e+00;
static const double  C4 = -2.549732539343734e+00;
static const double  C5 =  4.374664141464968e+00;
static const double  C6 =  2.938163982698783e+00;

static const double  D1 =  7.784695709041462e-03;
static const double  D2 =  3.224671290700398e-01;
static const double  D3 =  2.445134137142996e+00;
static const double  D4 =  3.754408661907416e+00;

static const double P_LOW =  0.02425;
static const double P_HIGH = 0.97575;  // P_high = 1 - p_low


long double ltqnorm(double p) {
  long double x = 0.0;  ///might be used uninitialized
  long double q, r, u, e;
  if ((0 < p )  && (p < P_LOW)){
    q = sqrt(-2.0*log(p));
    x = (((((C1*q+C2)*q+C3)*q+C4)*q+C5)*q+C6) / ((((D1*q+D2)*q+D3)*q+D4)*q+1);
  }
  else{
    if ((P_LOW <= p) && (p <= P_HIGH)){
      q = p - 0.5;
      r = q*q;
      x = (((((A1*r+A2)*r+A3)*r+A4)*r+A5)*r+A6)*q /(((((B1*r+B2)*r+B3)*r+B4)*r+B5)*r+1);
    }
    else{
      if ((P_HIGH < p)&&(p < 1)){
        q = sqrt(-2.0*log(1-p));
        x = -(((((C1*q+C2)*q+C3)*q+C4)*q+C5)*q+C6) / ((((D1*q+D2)*q+D3)*q+D4)*q+1);
      }
    }
  }

  //If you are compiling this under UNIX OR LINUX, you may uncomment this block for better accuracy.
  if(( 0 < p)&&(p < 1)) {
    e = 0.5 * erfc(-x/sqrt(2.0)) - p;
    u = e * sqrt(2.0*M_PI) * exp(x*x/2.0);
    x = x - u/(1 + x*u/2.0);
  }
    
  return x;
}


long double computeConfidenceInterval(double probability) {
  return ltqnorm((probability + 1) / 2);
}


