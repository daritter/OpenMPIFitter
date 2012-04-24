#include <iostream>
#include <fstream>
#include <sstream>

#include "tatami/tatami.h"

#include "Minuit2/MnUserParameters.h"

#include "TDirectory.h"

#include "cheb.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

// Chebyshev Polynomials
double chebyshev1( const double& x, const double& c )
{
  const double cheb_pdf =
    1.0 +
    (c*x);

  return cheb_pdf;
}

double norm_chebyshev1( const double& x_ll, const double& x_ul,
			const double& c )
{
  const double x_ll2 = x_ll*x_ll;

  const double x_ul2 = x_ul*x_ul;

  const double int_cheb_pdf =
    x_ul - x_ll +
    (c*0.5*x_ul2) - (c*0.5*x_ll2);

  return int_cheb_pdf;
}

double chebyshev2( const double& x, const std::vector<double>& c )
{
  if( c.size() != 2 )
    {
      std::cout << "ERROR: Number of coefficients incompatible with order. Check code." << std::endl;
      exit(1);
    }

  const double x2 = x*x;

  const double cheb_pdf =
    1.0 +
    (c.at(0)*x) +
    (c.at(1)*((2.0*x2) - 1.0));

  return cheb_pdf;
}

double norm_chebyshev2( const double& x_ll, const double& x_ul,
			const std::vector<double>& c )
{
  if( c.size() != 2 )
    {
      std::cout << "ERROR: Number of coefficients incompatible with order. Check code." << std::endl;
      exit(1);
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

double chebyshev3( const double& x, const std::vector<double>& c )
{
  if( c.size() != 3 )
    {
      std::cout << "ERROR: Number of coefficients incompatible with order. Check code." << std::endl;
      exit(1);
    }

  const double x2 = x*x;
  const double x3 = x2*x;

  const double cheb_pdf =
    1.0 +
    (c.at(0)*x) +
    (c.at(1)*((2.0*x2) - 1.0)) +
    (c.at(2)*((4.0*x3) - (3.0*x)));

  return cheb_pdf;
}

double norm_chebyshev3( const double& x_ll, const double& x_ul,
			const std::vector<double>& c )
{
  if( c.size() != 3 )
    {
      std::cout << "ERROR: Number of coefficients incompatible with order. Check code." << std::endl;
      exit(1);
    }

  const double x_ll2 = x_ll*x_ll;
  const double x_ll3 = x_ll2*x_ll;
  const double x_ll4 = x_ll3*x_ll;

  const double x_ul2 = x_ul*x_ul;
  const double x_ul3 = x_ul2*x_ul;
  const double x_ul4 = x_ul3*x_ul;

  const double int_cheb_pdf =
    x_ul - x_ll +
    (c.at(0)*0.5*x_ul2) - (c.at(0)*0.5*x_ll2) +
    (c.at(1)*((2.0/3.0*x_ul3) - x_ul)) - (c.at(1)*((2.0/3.0*x_ll3) - x_ll)) +
    (c.at(2)*(x_ul4 - (3.0/2.0*x_ul2))) - (c.at(2)*(x_ll4 - (3.0/2.0*x_ll2)));

  return int_cheb_pdf;
}

double chebyshev4( const double& x, const std::vector<double>& c )
{
  if( c.size() != 4 )
    {
      std::cout << "ERROR: Number of coefficients incompatible with order. Check code." << std::endl;
      exit(1);
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

double norm_chebyshev4( const double& x_ll, const double& x_ul,
			const std::vector<double>& c )
{
  if( c.size() != 4 )
    {
      std::cout << "ERROR: Number of coefficients incompatible with order. Check code." << std::endl;
      exit(1);
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

double chebyshev4( const double& xraw, const std::vector<double>& c,
		   const double& offset )
{
  if( c.size() != 4 )
    {
      std::cout << "ERROR: Number of coefficients incompatible with order. Check code." << std::endl;
      exit(1);
    }

  const double x  = xraw-offset;
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

double norm_chebyshev4( const double& xraw_ll, const double& xraw_ul,
			const std::vector<double>& c, const double& offset )
{
  if( c.size() != 4 )
    {
      std::cout << "ERROR: Number of coefficients incompatible with order. Check code." << std::endl;
      exit(1);
    }

  const double x_ll = xraw_ll - offset;
  const double x_ll2 = x_ll*x_ll;
  const double x_ll3 = x_ll2*x_ll;
  const double x_ll4 = x_ll3*x_ll;
  const double x_ll5 = x_ll4*x_ll;

  const double x_ul = xraw_ul - offset;
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

double chebyshev6( const double& x, const std::vector<double>& c )
{
  if( c.size() != 6 )
    {
      std::cout << "ERROR: Number of coefficients incompatible with order. Check code." << std::endl;
      exit(1);
    }

  const double x2 = x*x;
  const double x3 = x2*x;
  const double x4 = x3*x;
  const double x5 = x4*x;
  const double x6 = x5*x;

  const double cheb_pdf =
    1.0 +
    (c.at(0)*x) +
    (c.at(1)*((2.0*x2) - 1.0)) +
    (c.at(2)*((4.0*x3) - (3.0*x))) +
    (c.at(3)*((8.0*x4) - (8.0*x2) + 1.0)) +
    (c.at(4)*((16.0*x5) - (20.0*x3) + (5.0*x))) +
    (c.at(5)*((32.0*x6) - (48.0*x4) + (18.0*x2) - 1.0));

  return cheb_pdf;
}

double norm_chebyshev6( const double& x_ll, const double& x_ul,
			const std::vector<double>& c )
{
  if( c.size() != 6 )
    {
      std::cout << "ERROR: Number of coefficients incompatible with order. Check code." << std::endl;
      exit(1);
    }

  const double x_ll2 = x_ll*x_ll;
  const double x_ll3 = x_ll2*x_ll;
  const double x_ll4 = x_ll3*x_ll;
  const double x_ll5 = x_ll4*x_ll;
  const double x_ll6 = x_ll5*x_ll;
  const double x_ll7 = x_ll6*x_ll;

  const double x_ul2 = x_ul*x_ul;
  const double x_ul3 = x_ul2*x_ul;
  const double x_ul4 = x_ul3*x_ul;
  const double x_ul5 = x_ul4*x_ul;
  const double x_ul6 = x_ul5*x_ul;
  const double x_ul7 = x_ul6*x_ul;

  const double int_cheb_pdf =
    x_ul - x_ll +
    (c.at(0)*0.5*x_ul2) - (c.at(0)*0.5*x_ll2) +
    (c.at(1)*((2.0/3.0*x_ul3) - x_ul)) - (c.at(1)*((2.0/3.0*x_ll3) - x_ll)) +
    (c.at(2)*(x_ul4 - (3.0/2.0*x_ul2))) - (c.at(2)*(x_ll4 - (3.0/2.0*x_ll2))) +
    (c.at(3)*((8.0/5.0*x_ul5) - (8.0/3.0*x_ul3) + x_ul)) -
    (c.at(3)*((8.0/5.0*x_ll5) - (8.0/3.0*x_ll3) + x_ll)) +
    (c.at(4)*((8.0/3.0*x_ul6) - (5.0*x_ul4) + (5.0/2.0*x_ul2))) -
    (c.at(4)*((8.0/3.0*x_ll6) - (5.0*x_ll4) + (5.0/2.0*x_ll2))) +
    (c.at(5)*((32.0/7.0*x_ul7) - (48.0/5.0*x_ul5) + (6.0*x_ul3) - x_ul)) -
    (c.at(5)*((32.0/7.0*x_ll7) - (48.0/5.0*x_ll5) + (6.0*x_ll3) - x_ll));

  return int_cheb_pdf;
}

double chebyshev8( const double& x, const std::vector<double>& c )
{
  if( c.size() != 8 )
    {
      std::cout << "ERROR: Number of coefficients incompatible with order. Check code." << std::endl;
      exit(1);
    }

  const double x2 = x*x;
  const double x3 = x2*x;
  const double x4 = x3*x;
  const double x5 = x4*x;
  const double x6 = x5*x;
  const double x7 = x6*x;
  const double x8 = x7*x;

  const double cheb_pdf =
    1.0 +
    (c.at(0)*x) +
    (c.at(1)*((2.0*x2) - 1.0)) +
    (c.at(2)*((4.0*x3) - (3.0*x))) +
    (c.at(3)*((8.0*x4) - (8.0*x2) + 1.0)) +
    (c.at(4)*((16.0*x5) - (20.0*x3) + (5.0*x))) +
    (c.at(5)*((32.0*x6) - (48.0*x4) + (18.0*x2) - 1.0)) +
    (c.at(6)*((64.0*x7) - (112.0*x5) + (56.0*x3) - (7.0*x))) +
    (c.at(7)*((128.0*x8) - (256.0*x6) + (160.0*x4) - (32.0*x2) + 1.0));

  return cheb_pdf;
}

double norm_chebyshev8( const double& x_ll, const double& x_ul,
			const std::vector<double>& c )
{
  if( c.size() != 8 )
    {
      std::cout << "ERROR: Number of coefficients incompatible with order. Check code." << std::endl;
      exit(1);
    }

  const double x_ll2 = x_ll*x_ll;
  const double x_ll3 = x_ll2*x_ll;
  const double x_ll4 = x_ll3*x_ll;
  const double x_ll5 = x_ll4*x_ll;
  const double x_ll6 = x_ll5*x_ll;
  const double x_ll7 = x_ll6*x_ll;
  const double x_ll8 = x_ll7*x_ll;
  const double x_ll9 = x_ll8*x_ll;

  const double x_ul2 = x_ul*x_ul;
  const double x_ul3 = x_ul2*x_ul;
  const double x_ul4 = x_ul3*x_ul;
  const double x_ul5 = x_ul4*x_ul;
  const double x_ul6 = x_ul5*x_ul;
  const double x_ul7 = x_ul6*x_ul;
  const double x_ul8 = x_ul7*x_ul;
  const double x_ul9 = x_ul8*x_ul;

  const double int_cheb_pdf =
    x_ul - x_ll +
    (c.at(0)*0.5*x_ul2) - (c.at(0)*0.5*x_ll2) +
    (c.at(1)*((2.0/3.0*x_ul3) - x_ul)) - (c.at(1)*((2.0/3.0*x_ll3) - x_ll)) +
    (c.at(2)*(x_ul4 - (3.0/2.0*x_ul2))) - (c.at(2)*(x_ll4 - (3.0/2.0*x_ll2))) +
    (c.at(3)*((8.0/5.0*x_ul5) - (8.0/3.0*x_ul3) + x_ul)) -
    (c.at(3)*((8.0/5.0*x_ll5) - (8.0/3.0*x_ll3) + x_ll)) +
    (c.at(4)*((8.0/3.0*x_ul6) - (5.0*x_ul4) + (5.0/2.0*x_ul2))) -
    (c.at(4)*((8.0/3.0*x_ll6) - (5.0*x_ll4) + (5.0/2.0*x_ll2))) +
    (c.at(5)*((32.0/7.0*x_ul7) - (48.0/5.0*x_ul5) + (6.0*x_ul3) - x_ul)) -
    (c.at(5)*((32.0/7.0*x_ll7) - (48.0/5.0*x_ll5) + (6.0*x_ll3) - x_ll)) +
    (c.at(6)*((8.0*x_ul8) - (56.0/3.0*x_ul6) + (14.0*x_ul4) -
	      (7.0/2.0*x_ul2))) -
    (c.at(6)*((8.0*x_ll8) - (56.0/3.0*x_ll6) + (14.0*x_ll4) -
	      (7.0/2.0*x_ll2))) +
    (c.at(7)*((128.0/9.0*x_ul9) - (256.0/7.0*x_ul7) + (32.0*x_ul5) -
	      (32.0/3.0*x_ul3) + x_ul)) -
    (c.at(7)*((128.0/9.0*x_ll9) - (256.0/7.0*x_ll7) + (32.0*x_ll5) -
	      (32.0/3.0*x_ll3) + x_ll));

  return int_cheb_pdf;
}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
