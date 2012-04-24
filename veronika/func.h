#ifndef Func_H_
#define Func_H_

#include "belle.h"
#include "tatami/tatami.h"

#include "Data.h"

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
		   
double pipiDtQ( const Data& data, const double& dt, const int& q, const std::vector<double>& par, unsigned int par0=90 );
		  
double pipiDtQ_kolja( const Data& data, const double& dt, const int& q, const std::vector<double>& par, unsigned int par0=0 );

double pipiDtQ_ver( const Data& data, const double& dt, const int& q, const std::vector<double>& par, unsigned int par0=0 );

double pipiDtQ_ver_charged( const Data& data, const double& dt, const int& q, const std::vector<double>& par, unsigned int par0=0 );

void AddOutlier( const Data& data, const double& dt, double& pdf );

double QQbarDtQC( const Data& data, const double& dt, const int& q, /*const int& c,*/ const std::vector<double>& par, unsigned int par0 );

void AddOutDtQC( const Data& data, const double& dt, double& pdf );

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

#endif //Func_H_
