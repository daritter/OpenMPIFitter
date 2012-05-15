#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cfloat>
#include <cstdlib>

//#include "wtag.h"
#include "func.h"
//#include "constant.h"
#include "tatami/tatami.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif



double bigauss( const double& x, const double& mean,
		const double& sigmal, const double& sigmar )
{
  const double arg = x - mean;

  double coef;
  if( arg < 0.0 )
    coef = -0.5/sigmal/sigmal;
  else
    coef = -0.5/sigmar/sigmar;

  const double pdf = exp(coef*arg*arg);

  return pdf;
}

double norm_bigauss( const double& x_ll, const double& x_ul, const double& mean,
		     const double& sigmal, const double& sigmar )
{
  const double absSigmal = fabs(sigmal);
  const double absSigmar = fabs(sigmar);

  const double xscalel = sqrt(2.0)*absSigmal;
  const double xscaler = sqrt(2.0)*absSigmar;

  double int_pdf;
  if( x_ul < mean )
    int_pdf =
      absSigmal * ( erf((x_ul - mean)/xscalel) - erf((x_ll - mean)/xscalel) );
  else if( x_ll > mean )
    int_pdf =
      absSigmar * ( erf((x_ul - mean)/xscaler) - erf((x_ll - mean)/xscaler) );
  else
    int_pdf =
      absSigmar * erf((x_ul - mean)/xscaler) -
      absSigmal * erf((x_ll - mean)/xscalel);

  int_pdf *= sqrt(M_PI/2.0);

  return int_pdf;
}


// Crystal Ball
double crystalball( const double& x, const double& mean, const double& sigma,
		    const double& n, const double& alpha )
{
  double t = (x-mean)/sigma;
  if( alpha < 0.0 )
    t = -t;

  const double absAlpha = fabs(alpha);

  double pdf;
  if( t >= -absAlpha )
    pdf = exp(-0.5*t*t);
  else
    {
      const double A = pow(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
      const double B = (n/absAlpha) - absAlpha;

      pdf = A/pow(B - t, n);
    }
  return pdf;
}

double norm_crystalball( const double& x_ll, const double& x_ul,
			 const double& mean, const double& sigma,
			 const double& n, const double& alpha )
{
  const double absSigma = fabs(sigma);

  double t_ll = (x_ll-mean)/absSigma;
  double t_ul = (x_ul-mean)/absSigma;

  if( alpha < 0.0 )
    {
      const double tmp = t_ll;
      t_ll = -t_ul;
      t_ul = -tmp;
    }

  const double absAlpha = fabs(alpha);

  double int_pdf = 0.0;
  if( t_ll >= -absAlpha )
    int_pdf =
      absSigma*sqrt(M_PI/2.0)*( erf(t_ul/sqrt(2.0)) - erf(t_ll/sqrt(2.0)) );
  else if( t_ul <= -absAlpha )
    {
      const double A = pow(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
      const double B = (n/absAlpha) - absAlpha;

      int_pdf =
	A*absSigma/(1.0-n)*( (1.0/pow(B-t_ll, n-1.0)) -
			     (1.0/pow(B-t_ul, n-1.0)) );
    }
  else
    {
      const double A = pow(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
      const double B = (n/absAlpha) - absAlpha;

      const double term1 =
	A*absSigma/(1.0-n)*( (1.0/pow(B-t_ll,n-1.0)) -
			     (1.0/pow(n/absAlpha,n-1.0)) );

      const double term2 =
	absSigma*sqrt(M_PI/2.0)*( erf(t_ul/sqrt(2.0)) -
				  erf(-absAlpha/sqrt(2.0)) );

      int_pdf = term1 + term2;
    }
  return int_pdf;
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

double cheb2( const double& x, const std::vector<double>& c )
{
  if( c.size() != 2 )
    {
      std::cout << "ERROR: Number of coefficients incompatible with order. Check code." << std::endl;
      std::exit(1);
    }

  const double x2 = x*x;

  const double cheb_pdf =
    1.0 +
    (c.at(0)*x) +
    (c.at(1)*((2.0*x2) - 1.0));

  return cheb_pdf;
}

double norm_cheb2( const double& x_ll, const double& x_ul,
		   const std::vector<double>& c )
{
  if( c.size() != 2 )
    {
      std::cout << "ERROR: Number of coefficients incompatible with order. Check code." << std::endl;
      std::exit(1);
    }

  const double x_ll2 = x_ll*x_ll;
  const double x_ll3 = x_ll2*x_ll;

  const double x_ul2 = x_ul*x_ul;
  const double x_ul3 = x_ul2*x_ul;

  const double int_cheb_pdf =
    x_ul - x_ll +
    (c.at(0)*0.5*x_ul2) - (c.at(0)*0.5*x_ll2) +
    (c.at(1)*((2.0/3.0*x_ul3) - x_ul)) - (c.at(1)*((2.0/3.0*x_ll3) - x_ll));

  return int_cheb_pdf;
}

double cheb4( const double& x, const std::vector<double>& c )
{
  if( c.size() != 4 )
    {
      std::cout << "ERROR: Number of coefficients incompatible with order. Check code." << std::endl;
      std::exit(1);
    }

  const double x2 = x*x;
  const double x3 = x2*x;
  const double x4 = x3*x;

  const double cheb_pdf =
    1.0 +
    (c.at(0)*x) +
    (c.at(1)*((2.0*x2) - 1.0)) +
    (c.at(2)*((4.0*x3) - (3.0*x))) +
    (c.at(3)*((8.0*x4) - (8.0*x2) + 1.0));

  return cheb_pdf;
}

double norm_cheb4( const double& x_ll, const double& x_ul,
		   const std::vector<double>& c )
{
  if( c.size() != 4 )
    {
      std::cout << "ERROR: Number of coefficients incompatible with order. Check code." << std::endl;
      std::exit(1);
    }

  const double x_ll2 = x_ll*x_ll;
  const double x_ll3 = x_ll2*x_ll;
  const double x_ll4 = x_ll3*x_ll;
  const double x_ll5 = x_ll4*x_ll;

  const double x_ul2 = x_ul*x_ul;
  const double x_ul3 = x_ul2*x_ul;
  const double x_ul4 = x_ul3*x_ul;
  const double x_ul5 = x_ul4*x_ul;

  const double int_cheb_pdf =
    x_ul - x_ll +
    (c.at(0)*0.5*x_ul2) - (c.at(0)*0.5*x_ll2) +
    (c.at(1)*((2.0/3.0*x_ul3) - x_ul)) - (c.at(1)*((2.0/3.0*x_ll3) - x_ll)) +
    (c.at(2)*(x_ul4 - (3.0/2.0*x_ul2))) - (c.at(2)*(x_ll4 - (3.0/2.0*x_ll2))) +
    (c.at(3)*((8.0/5.0*x_ul5) - (8.0/3.0*x_ul3) + x_ul)) -
    (c.at(3)*((8.0/5.0*x_ll5) - (8.0/3.0*x_ll3) + x_ll));

  return int_cheb_pdf;
}


#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
