#ifndef Tools_H_
#define Tools_H_

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "Minuit2/FunctionMinimum.h"

#include "constant.h"
#include "Data.h"

using namespace ROOT::Minuit2;

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

// Load fit parameters
void load_par( std::ifstream& datafile, MnUserParameters& mn_param );

void compare_par( const MnUserParameters& par_in,
		  const FunctionMinimum& par_aus );

// phi2
double get_phi2eff( const double& Ccp, const double& DCcp,
		    const double& Scp, const double& DScp );
double get_phi2eff_err( const double& Ccp,  const double& DCcp,
			const double& Scp,  const double& DScp,
			const double& dCcp, const double& dDCcp,
			const double& dScp, const double& dDScp );

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

#endif //Tools_H_
