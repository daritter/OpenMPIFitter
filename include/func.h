#ifndef Func_H_
#define Func_H_

#include <vector>
#include "belle.h"
#include "tatami/tatami.h"

//#include "Data.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif


double bigauss( const double& x, const double& mean,
		const double& sigmal, const double& sigmar );
double norm_bigauss( const double& x_ll, const double& x_ul, const double& mean,
		     const double& sigmal, const double& sigmar );

//CRYSTAL BALL
double crystalball( const double& x, const double& mean, const double& sigma,
		    const double& n, const double& alpha );
double norm_crystalball( const double& x_ll, const double& x_ul,
			 const double& mean, const double& sigma,
			 const double& n, const double& alpha );

// ARGUS
double argus( const double& x,
	      const double& benergy, const double& a );
double norm_argus( const double& x_ll, const double& x_ul,
		   const double& benergy, const double& a );

// Chebyshev Polynomials
double cheb1( const double& x, const double& c );
double norm_cheb1( const double& x_ll, const double& x_ul,
		   const double& c );

double cheb2( const double& x, const std::vector<double>& c );
double norm_cheb2( const double& x_ll, const double& x_ul, const std::vector<double>& c );

double cheb4( const double& x, const std::vector<double>& c );
double norm_cheb4( const double& x_ll, const double& x_ul, const std::vector<double>& c );

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

#endif //Func_H_
