#ifndef Wtag_H_
#define Wtag_H_

#include "belle.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

//Wrong tag fraction
const unsigned int bins_rbin = 7;
const double rbin_bound[] = {0.0, 0.1, 0.25, 0.5, 0.625, 0.75, 0.875, 1.0};

const int set_rbin(const double r);

double set_r( const unsigned int& r_bin );

double set_wtag( const unsigned int& exp_no, const int& rbin,
		 const unsigned int& mc_type );
double set_dwtag( const unsigned int& exp_no, const int& rbin,
		  const unsigned int& mc_type );
double set_wtag_err( const unsigned int& exp_no, const int& rbin );
double set_dwtag_err( const unsigned int& exp_no, const int& rbin );

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

#endif //Wtag_H_
