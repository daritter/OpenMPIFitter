#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cfloat>

#include "func.h"

// Gaussian
double gaussian( const double& x, const double& m, const double& s )
{
  if(!finite(s))
    return 0;
  if(s==0.0)
    return FLT_MAX;
  double inv_s = 1.0/fabs(s);
  if(!finite(inv_s))
    return FLT_MAX;
  const double inv_sqrt_2pi =
    0.398942280401432702863218082711682654917240142822265625;
  const double dx = x-m;
  return inv_sqrt_2pi*inv_s*exp(-0.5*dx*dx*inv_s*inv_s);
}

double norm_gaussian( const double& x_ll, const double& x_ul,
		      const double& m, const double& s )
{
  double a_s = fabs(s);
  if(s==0.0) return 1;
  const double inv_sqrt2 = 0.707106781;
  double inv_s = 1.0/a_s;
  double x1 = (-x_ll+m)*inv_sqrt2*inv_s;
  double x2 = (x_ul-m)*inv_sqrt2*inv_s;
  return 1.0-0.5*erfc(x1)-0.5*erfc(x2);
}

// ARGUS
double argus( const double& x,
	      const double& benergy, const double& a )
{
  const double offset = 0.0;

  if ( x > benergy )
    return 0.0;

  double dz1 = x + offset;
  double dz2 = 1.0 - ( dz1/benergy ) * ( dz1/benergy );
  double dz3 = a * dz2;
  double argus = dz1 * sqrt(dz2) * exp(dz3);

  return argus;
}

double norm_argus( const double& x_ll, const double& x_ul,
		   const double& benergy, const double& a )
{
  const double offset = 0.0;

  double e2 = benergy * benergy;
  double ll2 = 1.0 - ( x_ll+offset ) * ( x_ll+offset ) / e2;
  double ul2 = 1.0 - ( x_ul+offset ) * ( x_ul+offset ) / e2;

  if( ll2 < 0.0 )
    return 0.0;
  if( ul2 < 0.0 )
    ul2 = 0.0;

  return 0.5*e2/a*((sqrt(ll2)*exp(a*ll2) - sqrt(ul2)*exp(a*ul2))
                   +0.5*sqrt(M_PI/fabs(a))
                   *(erfc(sqrt(fabs(a)*ll2)) - erfc(sqrt(fabs(a)*ul2))));
}

// Chebyshev Polynomials
double cheb1( const double& x, const double& c )
{
  const double cheb_pdf =
    1.0 +
    (c*x);

  return cheb_pdf;
}

double norm_cheb1( const double& x_ll, const double& x_ul,
		   const double& c )
{
  const double x_ll2 = x_ll*x_ll;

  const double x_ul2 = x_ul*x_ul;

  const double int_cheb_pdf =
    x_ul - x_ll +
    (c*0.5*x_ul2) - (c*0.5*x_ll2);

  return int_cheb_pdf;
}
