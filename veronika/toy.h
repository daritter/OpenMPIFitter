#ifndef Toy_H_
#define Toy_H_

#include "TRandom3.h"
#include "TH1.h"

#include "Data.h"
#include "tools.h"

using namespace ROOT::Minuit2;

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

// Random distribution generators
void rand_dist2D( TRandom3* rand,
		  double& x, double& h,
		  const double& xll, const double& xul,
		  const double& max );

void rand_dist3D( TRandom3* rand,
		  double& x, double& y, double& h,
		  const double& xll, const double& xul,
		  const double& yll, const double& yul,
		  const double& max );

void rand_dist4D( TRandom3* rand,
		  double& x, double& y, double& z, double& h,
		  const double& xll, const double& xul,
		  const double& yll, const double& yul,
		  const double& zll, const double& zul,
		  const double& max );
void rand_dist5D( TRandom3* rand,
		  double &a, double &b, double &c,
		  double &d, double &h,
		  const double& all, const double& aul,
		  const double& bll, const double& bul,
		  const double& cll, const double& cul,
		  const double& dll, const double& dul,
		  const double& max );
void rand_dist6D( TRandom3* rand,
		  double &a, double &b, double &c,
		  double &d, double &e, double &h,
		  const double& all, const double& aul,
		  const double& bll, const double& bul,
		  const double& cll, const double& cul,
		  const double& dll, const double& dul,
		  const double& ell, const double& eul,
		  const double& max );
unsigned int vector_rand( TRandom3* rand,
			  std::vector<unsigned int>& vec );
double vector_rand( TRandom3* rand,
		    std::vector<double>& vec );

Data vector_rand( TRandom3* rand, std::vector<Data>& vec );

double gen_hist1d( TRandom3* rand, const TH1D* h_ref );

double gen_r( TRandom3* rand, const TH1D* h_rbin );

double gen_costhetaB_sig( TRandom3* rand );
double gen_costhetaB_bkg( TRandom3* rand );

int gen_q( TRandom3* rand );
int gen_c( TRandom3* rand );

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

#endif //Toy_H_
