#ifndef Cheb_H_
#define Cheb_H_

using namespace ROOT::Minuit2;

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

// Chebyshev Polynomials
double chebyshev1( const double& x, const double& c );
double norm_chebyshev1( const double& x_ll, const double& x_ul,
			const double& c );
double chebyshev2( const double& x, const std::vector<double>& c );
double norm_chebyshev2( const double& x_ll, const double& x_ul,
			const std::vector<double>& c );
double chebyshev3( const double& x, const std::vector<double>& c );
double norm_chebyshev3( const double& x_ll, const double& x_ul,
			const std::vector<double>& c );
double chebyshev4( const double& x, const std::vector<double>& c );
double norm_chebyshev4( const double& x_ll, const double& x_ul,
			const std::vector<double>& c );
double chebyshev4( const double& xraw, const std::vector<double>& c,
		   const double& offset );
double norm_chebyshev4( const double& xraw_ll, const double& xraw_ul,
			const std::vector<double>& c, const double& offset );
double chebyshev6( const double& x, const std::vector<double>& c );
double norm_chebyshev6( const double& x_ll, const double& x_ul,
			const std::vector<double>& c );
double chebyshev8( const double& x, const std::vector<double>& c );
double norm_chebyshev8( const double& x_ll, const double& x_ul,
			const std::vector<double>& c );

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

#endif //Cheb_H_
