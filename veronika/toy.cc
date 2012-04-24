#include <iostream>
#include <fstream>
#include <sstream>

#include "tatami/tatami.h"

#include "Minuit2/MnUserParameters.h"

#include "TDirectory.h"

#include "hist.h"
#include "wtag.h"
#include "toy.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

// Random distribution generators
void rand_dist2D( TRandom3* rand,
		  double &x, double &h,
		  const double& xll, const double& xul,
		  const double& max )
{
  x = rand->Uniform(xll, xul);
  h = rand->Uniform(0.0, max);
}

void rand_dist3D( TRandom3* rand,
		  double &x, double &y, double &h,
		  const double& xll, const double& xul,
		  const double& yll, const double& yul,
		  const double& max )
{
  x = rand->Uniform(xll, xul);
  y = rand->Uniform(yll, yul);
  h = rand->Uniform(0.0, max);
}

void rand_dist4D( TRandom3* rand,
		  double &x, double &y, double &z, double &h,
		  const double& xll, const double& xul,
		  const double& yll, const double& yul,
		  const double& zll, const double& zul,
		  const double& max )
{
  x = rand->Uniform(xll, xul);
  y = rand->Uniform(yll, yul);
  z = rand->Uniform(zll, zul);
  h = rand->Uniform(0.0, max);
}

void rand_dist5D( TRandom3* rand,
		  double &a, double &b, double &c,
		  double &d, double &h,
		  const double& all, const double& aul,
		  const double& bll, const double& bul,
		  const double& cll, const double& cul,
		  const double& dll, const double& dul,
		  const double& max )
{
  a = rand->Uniform(all, aul);
  b = rand->Uniform(bll, bul);
  c = rand->Uniform(cll, cul);
  d = rand->Uniform(dll, dul);
  h = rand->Uniform(0.0, max);
}

void rand_dist6D( TRandom3* rand,
		  double &a, double &b, double &c,
		  double &d, double &e, double &h,
		  const double& all, const double& aul,
		  const double& bll, const double& bul,
		  const double& cll, const double& cul,
		  const double& dll, const double& dul,
		  const double& ell, const double& eul,
		  const double& max )
{
  a = rand->Uniform(all, aul);
  b = rand->Uniform(bll, bul);
  c = rand->Uniform(cll, cul);
  d = rand->Uniform(dll, dul);
  e = rand->Uniform(ell, eul);
  h = rand->Uniform(0.0, max);
}

unsigned int vector_rand( TRandom3* rand,
			  std::vector<unsigned int>& vec )
{
  //Pick random entry
  const unsigned int entry =
    static_cast<unsigned int>(round( rand->Rndm()*
				     static_cast<double>(vec.size()-1) ));

  const unsigned int data = vec.at(entry);
  vec.erase( vec.begin()+entry );

  return data;
}

double vector_rand( TRandom3* rand,
		    std::vector<double>& vec )
{
  //Pick random entry
  const unsigned int entry =
    static_cast<unsigned int>(round( rand->Rndm()*
				     static_cast<double>(vec.size()-1) ));

  const double data = vec.at(entry);
  vec.erase( vec.begin()+entry );

  return data;
}

Data vector_rand( TRandom3* rand, std::vector<Data>& vec )
{
  //Pick random entry
  const unsigned int entry =
    static_cast<unsigned int>(round( rand->Rndm()*
				     static_cast<double>(vec.size()-1) ));

  const Data data = vec.at(entry);
  vec.erase( vec.begin()+entry );

  return data;
}

double gen_hist1d( TRandom3* rand, const TH1D* h_ref )
{
  const double height = h_ref->GetMaximum();

  for(;;)
    {
      double x, h;

      rand_dist2D( rand, x, h,
		   h_ref->GetXaxis()->GetXmin(), h_ref->GetXaxis()->GetXmax(),
		   height );

      const double pdf = hist1d( x, h_ref );
      if( h < pdf )
	return x;
    }
}

double gen_r( TRandom3* rand, const TH1D* h_rbin )
{
  //Account for variable bin width
  TH1D h_rbin_( "h_rbin_", "", bins_rbin, rbin_bound );
  h_rbin_.Add(h_rbin);

  for( unsigned int i=0; i<bins_rbin; i++ )
    h_rbin_.
      SetBinContent( i+1,
		     h_rbin_.GetBinContent(i+1)*h_rbin_.GetBinWidth(1)/
		     h_rbin_.GetBinWidth(i+1) );

  const double height_max = 1.1*h_rbin_.GetMaximum();
  for(;;)
    {
      double x, h;

      rand_dist2D( rand, x, h,
		   0.0, 1.0,
		   height_max );

      const double pdf = hist1d( x, &h_rbin_ );
      if( h < pdf )
	return x;
    }
}

double gen_costhetaB_sig( TRandom3* rand )
{
  const double height_max = 1.0;

  double costhetaB_gen, height_gen;
  for(;;)
    {
      rand_dist2D( rand, costhetaB_gen, height_gen,
		   -1.0, 1.0, height_max );
      if( height_gen < 1.0-(costhetaB_gen*costhetaB_gen) )
	break;
    }
  return costhetaB_gen;
}

double gen_costhetaB_bkg( TRandom3* rand )
{
  return rand->Uniform(-1.0, 1.0);
}

int gen_q( TRandom3* rand )
{
  int q;
  if( rand->Rndm() < 0.5 )
    q = -1;
  else
    q = 1;

  return q;
}

int gen_c( TRandom3* rand )
{
  return gen_q(rand);
}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
