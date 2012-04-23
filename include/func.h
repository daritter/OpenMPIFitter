#ifndef Func_H_
#define Func_H_

// Gaussian
double gaussian( const double& x, const double& m, const double& s );
double norm_gaussian( const double& x_ll, const double& x_ul,
		      const double& m, const double& s );

// ARGUS
double argus( const double& x,
	      const double& benergy, const double& a );
double norm_argus( const double& x_ll, const double& x_ul,
		   const double& benergy, const double& a );

// Chebyshev Polynomials
double cheb1( const double& x, const double& c );
double norm_cheb1( const double& x_ll, const double& x_ul,
		   const double& c );

#endif //Func_H_
