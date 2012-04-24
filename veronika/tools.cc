#include <iostream>
#include <fstream>
#include <sstream>

#include "tatami/tatami.h"

#include "TDirectory.h"
#include "Minuit2/MnUserParameters.h"

#include "tools.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

// Load fit parameters
void load_par( std::ifstream& datafile, MnUserParameters& mn_param )
{
  if( datafile.is_open() == 0 )
    {
      std::cerr << "ERROR: Datafile could not be opened."
                << std::endl;
      exit(1);
    }
  std::string dummy_name;
  double dummy_par, dummy_err;
  datafile >> dummy_name >> dummy_par >> dummy_err;
  while( datafile.eof() == 0 )
    {
      mn_param.Add(dummy_name.c_str(), dummy_par, dummy_err);
      datafile >> dummy_name >> dummy_par >> dummy_err;
    }
  datafile.close();
}

void compare_par( const MnUserParameters& par_in,
		  const FunctionMinimum& par_aus )
{
  for( unsigned int i=0; i<par_in.Params().size(); i++ )
    std::cout << par_in.Name(i) << ": "
	      << fabs(par_in.Params().at(i) - par_aus.UserParameters().Params().at(i))/par_aus.UserParameters().Errors().at(i) << "s" <<std::endl;
}

// phi2
double get_phi2eff( const double& Ccp, const double& DCcp,
		    const double& Scp, const double& DScp )
{
  const double phi2 =
    180.0/M_PI/4.0*
    (asin((Scp + DScp)/sqrt(1.0-((Ccp+DCcp)*(Ccp+DCcp)))) +
     asin((Scp - DScp)/sqrt(1.0-((Ccp-DCcp)*(Ccp-DCcp)))));

  std::cout << 180.0/M_PI/4.0*
    (asin((Scp + DScp)/sqrt(1.0-((Ccp+DCcp)*(Ccp+DCcp)))) +
     asin((Scp - DScp)/sqrt(1.0-((Ccp-DCcp)*(Ccp-DCcp)))))
	    << std::endl;
  std::cout << 180.0/M_PI/4.0*
    ((M_PI-asin((Scp + DScp)/sqrt(1.0-((Ccp+DCcp)*(Ccp+DCcp))))) +
     asin((Scp - DScp)/sqrt(1.0-((Ccp-DCcp)*(Ccp-DCcp)))))
	    << std::endl;
  std::cout << 180.0/M_PI/4.0*
    (asin((Scp + DScp)/sqrt(1.0-((Ccp+DCcp)*(Ccp+DCcp)))) +
     (M_PI-asin((Scp - DScp)/sqrt(1.0-((Ccp-DCcp)*(Ccp-DCcp))))))
	    << std::endl;
  std::cout << 180.0/M_PI/4.0*
    ((M_PI-asin((Scp + DScp)/sqrt(1.0-((Ccp+DCcp)*(Ccp+DCcp))))) +
     (M_PI-asin((Scp - DScp)/sqrt(1.0-((Ccp-DCcp)*(Ccp-DCcp))))))
	    << std::endl;

  return phi2;
}

double get_phi2eff_err( const double& Ccp,  const double& DCcp,
			const double& Scp,  const double& DScp,
			const double& dCcp, const double& dDCcp,
			const double& dScp, const double& dDScp )
{
  const double g = 1.0-((Ccp+DCcp)*(Ccp+DCcp));
  const double h = 1.0-((Ccp-DCcp)*(Ccp-DCcp));

  const double x = (Scp + DScp)/sqrt(g);
  const double y = (Scp - DScp)/sqrt(h);
  const double x2 = x*x;
  const double y2 = y*y;

  const double dfdScp =
    (1.0/4.0/sqrt(1.0-x2)/sqrt(g)) +
    (1.0/4.0/sqrt(1.0-y2)/sqrt(h));
  const double dfdDScp =
    (1.0/4.0/sqrt(1.0-x2)/sqrt(g)) -
    (1.0/4.0/sqrt(1.0-y2)/sqrt(h));

  const double dfdCcp =
    (1.0/4.0/sqrt(1.0-x2)*(Scp+DScp)*(Ccp+DCcp)/sqrt(g)/g) +
    (1.0/4.0/sqrt(1.0-y2)*(Scp-DScp)*(Ccp-DCcp)/sqrt(h)/h);

  const double dfdDCcp =
    (1.0/4.0/sqrt(1.0-x2)*(Scp+DScp)*(Ccp+DCcp)/sqrt(g)/g) -
    (1.0/4.0/sqrt(1.0-y2)*(Scp-DScp)*(Ccp-DCcp)/sqrt(h)/h);

  std::cout << sqrt(g) << std::endl;
  std::cout << sqrt(h) << std::endl;
  std::cout << x2 << std::endl;
  std::cout << sqrt(1.0-x2) << std::endl;
  std::cout << sqrt(1.0-y2) << std::endl;
  std::cout << dfdScp << std::endl;
  std::cout << dfdDScp << std::endl;
  std::cout << dfdCcp << std::endl;
  std::cout << dfdDCcp << std::endl;

  const double dphi2eff =
    sqrt((dfdScp*dfdScp*dScp*dScp) +
	 (dfdDScp*dfdDScp*dDScp*dDScp) +
	 (dfdCcp*dfdCcp*dCcp*dCcp) +
	 (dfdDCcp*dfdDCcp*dDCcp*dDCcp));

  return 180.0/M_PI*dphi2eff;
}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
