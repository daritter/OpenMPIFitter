#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cfloat>

#include "wtag.h"
#include "func.h"
#include "constant.h"
#include "tatami/tatami.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif



double bigauss( const double& x, const double& mean,
		const double& sigmal, const double& sigmar )
{
  const double arg = x - mean;

  double coef;
  if( arg < 0.0 )
    coef = -0.5/sigmal/sigmal;
  else
    coef = -0.5/sigmar/sigmar;

  const double pdf = exp(coef*arg*arg);

  return pdf;
}

double norm_bigauss( const double& x_ll, const double& x_ul, const double& mean,
		     const double& sigmal, const double& sigmar )
{
  const double absSigmal = fabs(sigmal);
  const double absSigmar = fabs(sigmar);

  const double xscalel = sqrt(2.0)*absSigmal;
  const double xscaler = sqrt(2.0)*absSigmar;

  double int_pdf;
  if( x_ul < mean )
    int_pdf =
      absSigmal * ( erf((x_ul - mean)/xscalel) - erf((x_ll - mean)/xscalel) );
  else if( x_ll > mean )
    int_pdf =
      absSigmar * ( erf((x_ul - mean)/xscaler) - erf((x_ll - mean)/xscaler) );
  else
    int_pdf =
      absSigmar * erf((x_ul - mean)/xscaler) -
      absSigmal * erf((x_ll - mean)/xscalel);

  int_pdf *= sqrt(M_PI/2.0);

  return int_pdf;
}


// Crystal Ball
double crystalball( const double& x, const double& mean, const double& sigma,
		    const double& n, const double& alpha )
{
  double t = (x-mean)/sigma;
  if( alpha < 0.0 )
    t = -t;

  const double absAlpha = fabs(alpha);

  double pdf;
  if( t >= -absAlpha )
    pdf = exp(-0.5*t*t);
  else
    {
      const double A = pow(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
      const double B = (n/absAlpha) - absAlpha; 

      pdf = A/pow(B - t, n);
    }
  return pdf;
}

double norm_crystalball( const double& x_ll, const double& x_ul,
			 const double& mean, const double& sigma,
			 const double& n, const double& alpha )
{
  const double absSigma = fabs(sigma);

  double t_ll = (x_ll-mean)/absSigma;
  double t_ul = (x_ul-mean)/absSigma;

  if( alpha < 0.0 )
    {
      const double tmp = t_ll;
      t_ll = -t_ul;
      t_ul = -tmp;
    }

  const double absAlpha = fabs(alpha);

  double int_pdf = 0.0;
  if( t_ll >= -absAlpha )
    int_pdf =
      absSigma*sqrt(M_PI/2.0)*( erf(t_ul/sqrt(2.0)) - erf(t_ll/sqrt(2.0)) );
  else if( t_ul <= -absAlpha )
    {
      const double A = pow(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
      const double B = (n/absAlpha) - absAlpha; 

      int_pdf =
	A*absSigma/(1.0-n)*( (1.0/pow(B-t_ll, n-1.0)) -
			     (1.0/pow(B-t_ul, n-1.0)) );
    }
  else
    {
      const double A = pow(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
      const double B = (n/absAlpha) - absAlpha; 

      const double term1 =
	A*absSigma/(1.0-n)*( (1.0/pow(B-t_ll,n-1.0)) -
			     (1.0/pow(n/absAlpha,n-1.0)) );

      const double term2 =
	absSigma*sqrt(M_PI/2.0)*( erf(t_ul/sqrt(2.0)) -
				  erf(-absAlpha/sqrt(2.0)) );

      int_pdf = term1 + term2;
    }
  return int_pdf;
}

// ARGUS
double argus( const double& x,
	      const double& benergy, const double& a )
{
  const double offset = 0.0;

  if ( x > benergy )
    return 0.0;

  double dz1 = x + offset;
  double dz2 = 1.0 - ( dz1/benergy ) * ( dz1/benergy );
  double dz3 = a * dz2;
  double argus = dz1 * sqrt(dz2) * exp(dz3);

  return argus;
}

double norm_argus( const double& x_ll, const double& x_ul,
		   const double& benergy, const double& a )
{
  const double offset = 0.0;

  double e2 = benergy * benergy;
  double ll2 = 1.0 - ( x_ll+offset ) * ( x_ll+offset ) / e2;
  double ul2 = 1.0 - ( x_ul+offset ) * ( x_ul+offset ) / e2;

  if( ll2 < 0.0 )
    return 0.0;
  if( ul2 < 0.0 )
    ul2 = 0.0;

  return 0.5*e2/a*((sqrt(ll2)*exp(a*ll2) - sqrt(ul2)*exp(a*ul2))
                   +0.5*sqrt(M_PI/fabs(a))
                   *(erfc(sqrt(fabs(a)*ll2)) - erfc(sqrt(fabs(a)*ul2))));
}

// Chebyshev Polynomials
double cheb1( const double& x, const double& c )
{
  const double cheb_pdf =
    1.0 +
    (c*x);

  return cheb_pdf;
}

double norm_cheb1( const double& x_ll, const double& x_ul,
		   const double& c )
{
  const double x_ll2 = x_ll*x_ll;

  const double x_ul2 = x_ul*x_ul;

  const double int_cheb_pdf =
    x_ul - x_ll +
    (c*0.5*x_ul2) - (c*0.5*x_ll2);

  return int_cheb_pdf;
}

double cheb2( const double& x, const std::vector<double>& c )
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

double norm_cheb2( const double& x_ll, const double& x_ul,
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

double cheb4( const double& x, const std::vector<double>& c )
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

double norm_cheb4( const double& x_ll, const double& x_ul,
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

double pipiDtQ( const Data& data,
			  const double& dt, const int& q,
			  const std::vector<double>& par,
			  unsigned int par0 )
{
 
  double Acp_pipi, Scp_pipi;

      Acp_pipi = par[par0+0];
      Scp_pipi = par[par0+1];

  const dtres_param_t* const dtres_param =
    get_dtres_param( data.GetExpNo(), data.GetMC() );
/*

  //Lifetime component
  double life_pdf;
    life_pdf =
      EfRkRdetRnp_fullrec( dt, data.GetBtype(), data.GetTau(),
			   data.GetAk(), data.GetCk(),
			   data.GetRecVNtrk(), data.GetRecVZerr(),
			   data.GetRecVChisq(), data.GetRecVNdf(),
			   data.GetTagVNtrk(), data.GetTagVZerr(),
			   data.GetTagVChisq(), data.GetTagVNdf(),
			   data.GetTagVIsl(), dtres_param );

  double int_life_pdf;
    int_life_pdf =
      norm_EfRkRdetRnp_fullrec( h_dt_ll, h_dt_ul, data.GetBtype(), 
				data.GetTau(), data.GetAk(), data.GetCk(),
				data.GetRecVNtrk(), data.GetRecVZerr(),
				data.GetRecVChisq(), data.GetRecVNdf(),
				data.GetTagVNtrk(), data.GetTagVZerr(),
				data.GetTagVChisq(), data.GetTagVNdf(),
				data.GetTagVIsl(), dtres_param );

  //Acp component
  double cos_pdf;
    cos_pdf =
      0.5 / data.GetTau() *
      MfRkRdetRnp_fullrec( dt, data.GetBtype(), data.GetTau(),
			   data.GetDeltam(), data.GetAk(), data.GetCk(),
			   data.GetRecVNtrk(), data.GetRecVZerr(),
			   data.GetRecVChisq(), data.GetRecVNdf(),
			   data.GetTagVNtrk(), data.GetTagVZerr(),
			   data.GetTagVChisq(), data.GetTagVNdf(),
			   data.GetTagVIsl(), dtres_param );

  //Scp component
  double sin_pdf;
      sin_pdf =
      0.5 / data.GetTau() *
      AfRkRdetRnp_fullrec( dt, data.GetBtype(), data.GetTau(),
			   data.GetDeltam(), data.GetAk(), data.GetCk(),
			   data.GetRecVNtrk(), data.GetRecVZerr(),
			   data.GetRecVChisq(), data.GetRecVNdf(),
			   data.GetTagVNtrk(), data.GetTagVZerr(),
			   data.GetTagVChisq(), data.GetTagVNdf(),
			   data.GetTagVIsl(), dtres_param );
			   
*/

  //Lifetime component
  const double life_pdf = data.GetEfRkRdetRnp();
  const double int_life_pdf = data.GetNormEfRkRdetRnp();

  //Acp component
  const double cos_pdf = data.GetMfRkRdetRnp();

  //Scp component
  const double sin_pdf = data.GetAfRkRdetRnp();

  //Wrong tag fraction
  const int r_bin = set_rbin( data.GetR() );
  const double wt_w = set_wtag( data.GetExpNo(), data.GetMC(), r_bin );
  const double wt_dw = set_dwtag( data.GetExpNo(), data.GetMC(), r_bin );

  //pipi PDF
  double pipi_pdf =
    (life_pdf * ( 1.0 - (q*wt_dw) )) + 
    (( q * (1.0 - (2.0*wt_w)) ) *
     ((Acp_pipi*cos_pdf)+(Scp_pipi*sin_pdf)));

  //Normalise over q
  pipi_pdf /= 2.0*int_life_pdf;

  AddOutlier(data, dt, pipi_pdf);

   return pipi_pdf;
}


double pipiDtQ_kolja( const Data& data,
			  const double& dt, const int& q,
			  const std::vector<double>& par,
			  unsigned int par0 )
{
 
  double Acp_pipi, Scp_pipi;

      Acp_pipi = par[par0+0];
      Scp_pipi = par[par0+1];

  const dtres_param_t* const dtres_param =
    get_dtres_param( data.GetExpNo(), data.GetMC() );

  //Lifetime component
  double life_pdf;
    life_pdf =
      EfRkRdetRnp_fullrec( dt, data.GetBtype(), data.GetTau(),
			   data.GetAk(), data.GetCk(),
			   data.GetRecVNtrk(), data.GetRecVZerr(),
			   data.GetRecVChisq(), data.GetRecVNdf(),
			   data.GetTagVNtrk(), data.GetTagVZerr(),
			   data.GetTagVChisq(), data.GetTagVNdf(),
			   data.GetTagVIsl(), dtres_param );

  double int_life_pdf;
    int_life_pdf =
      norm_EfRkRdetRnp_fullrec( h_dt_ll, h_dt_ul, data.GetBtype(), 
				data.GetTau(), data.GetAk(), data.GetCk(),
				data.GetRecVNtrk(), data.GetRecVZerr(),
				data.GetRecVChisq(), data.GetRecVNdf(),
				data.GetTagVNtrk(), data.GetTagVZerr(),
				data.GetTagVChisq(), data.GetTagVNdf(),
				data.GetTagVIsl(), dtres_param );

  //Acp component
  double cos_pdf;
    cos_pdf =
      0.5 / data.GetTau() *
      MfRkRdetRnp_fullrec( dt, data.GetBtype(), data.GetTau(),
			   data.GetDeltam(), data.GetAk(), data.GetCk(),
			   data.GetRecVNtrk(), data.GetRecVZerr(),
			   data.GetRecVChisq(), data.GetRecVNdf(),
			   data.GetTagVNtrk(), data.GetTagVZerr(),
			   data.GetTagVChisq(), data.GetTagVNdf(),
			   data.GetTagVIsl(), dtres_param );

  //Scp component
  double sin_pdf;
      sin_pdf =
      0.5 / data.GetTau() *
      AfRkRdetRnp_fullrec( dt, data.GetBtype(), data.GetTau(),
			   data.GetDeltam(), data.GetAk(), data.GetCk(),
			   data.GetRecVNtrk(), data.GetRecVZerr(),
			   data.GetRecVChisq(), data.GetRecVNdf(),
			   data.GetTagVNtrk(), data.GetTagVZerr(),
			   data.GetTagVChisq(), data.GetTagVNdf(),
			   data.GetTagVIsl(), dtres_param );
			   
  //Wrong tag fraction
  const int r_bin = set_rbin( data.GetR() );
  const double wt_w = set_wtag( data.GetExpNo(), data.GetMC(), r_bin );
  const double wt_dw = set_dwtag( data.GetExpNo(), data.GetMC(), r_bin );

  //pipi PDF
  double pipi_pdf =
    (life_pdf * ( 1.0 - (q*wt_dw) )) + 
    (( q * (1.0 - (2.0*wt_w)) ) *
     ((Acp_pipi*cos_pdf)+(Scp_pipi*sin_pdf)));

  //Normalise over q
  pipi_pdf /= 2.0*int_life_pdf;

  AddOutlier(data, dt, pipi_pdf);

   return pipi_pdf;
}

double pipiDtQ_ver( const Data& data,
			  const double& dt, const int& q,
			  const std::vector<double>& par,
			  unsigned int par0 )
{
 
  double Acp_pipi, Scp_pipi, Tau;

      Acp_pipi = par[par0+0];
      Scp_pipi = par[par0+1];
      Tau = par[par0+2];

  const dtres_param_t* const dtres_param =
    get_dtres_param( data.GetExpNo(), data.GetMC() );

  //Lifetime component
  double life_pdf;
    life_pdf =
      EfRkRdetRnp_fullrec( dt, data.GetBtype(), Tau,
			   data.GetAk(), data.GetCk(),
			   data.GetRecVNtrk(), data.GetRecVZerr(),
			   data.GetRecVChisq(), data.GetRecVNdf(),
			   data.GetTagVNtrk(), data.GetTagVZerr(),
			   data.GetTagVChisq(), data.GetTagVNdf(),
			   data.GetTagVIsl(), dtres_param );

  double int_life_pdf;
    int_life_pdf =
      norm_EfRkRdetRnp_fullrec( h_dt_ll, h_dt_ul, data.GetBtype(), 
				Tau, data.GetAk(), data.GetCk(),
				data.GetRecVNtrk(), data.GetRecVZerr(),
				data.GetRecVChisq(), data.GetRecVNdf(),
				data.GetTagVNtrk(), data.GetTagVZerr(),
				data.GetTagVChisq(), data.GetTagVNdf(),
				data.GetTagVIsl(), dtres_param );

  //Acp component
  double cos_pdf;
    cos_pdf =
      0.5 / Tau *
      MfRkRdetRnp_fullrec( dt, data.GetBtype(), Tau,
			   data.GetDeltam(), data.GetAk(), data.GetCk(),
			   data.GetRecVNtrk(), data.GetRecVZerr(),
			   data.GetRecVChisq(), data.GetRecVNdf(),
			   data.GetTagVNtrk(), data.GetTagVZerr(),
			   data.GetTagVChisq(), data.GetTagVNdf(),
			   data.GetTagVIsl(), dtres_param );

  //Scp component
  double sin_pdf;
      sin_pdf =
      0.5 / Tau *
      AfRkRdetRnp_fullrec( dt, data.GetBtype(), Tau,
			   data.GetDeltam(), data.GetAk(), data.GetCk(),
			   data.GetRecVNtrk(), data.GetRecVZerr(),
			   data.GetRecVChisq(), data.GetRecVNdf(),
			   data.GetTagVNtrk(), data.GetTagVZerr(),
			   data.GetTagVChisq(), data.GetTagVNdf(),
			   data.GetTagVIsl(), dtres_param );
			   
  //Wrong tag fraction
  const int r_bin = set_rbin( data.GetR() );
  const double wt_w = set_wtag( data.GetExpNo(), data.GetMC(), r_bin );
  const double wt_dw = set_dwtag( data.GetExpNo(), data.GetMC(), r_bin );

  //pipi PDF
  double pipi_pdf =
    (life_pdf * ( 1.0 - (q*wt_dw) )) + 
    (( q * (1.0 - (2.0*wt_w)) ) *
     ((Acp_pipi*cos_pdf)+(Scp_pipi*sin_pdf)));

  //Normalise over q
  pipi_pdf /= 2.0*int_life_pdf;

  AddOutlier(data, dt, pipi_pdf);

   return pipi_pdf;
}

double pipiDtQ_ver_charged( const Data& data,
			  const double& dt, const int& q,
			  const std::vector<double>& par,
			  unsigned int par0 )
{
 
  double Tau;

      Tau = par[par0+0];

  const dtres_param_t* const dtres_param =
    get_dtres_param( data.GetExpNo(), data.GetMC() );

  //Lifetime component
  double life_pdf;
    life_pdf =
      EfRkRdetRnp_fullrec( dt, 1., Tau,
			   data.GetAk(), data.GetCk(),
			   data.GetRecVNtrk(), data.GetRecVZerr(),
			   data.GetRecVChisq(), data.GetRecVNdf(),
			   data.GetTagVNtrk(), data.GetTagVZerr(),
			   data.GetTagVChisq(), data.GetTagVNdf(),
			   data.GetTagVIsl(), dtres_param );

  double int_life_pdf;
    int_life_pdf =
      norm_EfRkRdetRnp_fullrec( h_dt_ll, h_dt_ul, 1., 
				Tau, data.GetAk(), data.GetCk(),
				data.GetRecVNtrk(), data.GetRecVZerr(),
				data.GetRecVChisq(), data.GetRecVNdf(),
				data.GetTagVNtrk(), data.GetTagVZerr(),
				data.GetTagVChisq(), data.GetTagVNdf(),
				data.GetTagVIsl(), dtres_param );

  //Acp component
  double cos_pdf;
    cos_pdf =
      0.5 / Tau *
      MfRkRdetRnp_fullrec( dt, 1., Tau,
			   data.GetDeltam(), data.GetAk(), data.GetCk(),
			   data.GetRecVNtrk(), data.GetRecVZerr(),
			   data.GetRecVChisq(), data.GetRecVNdf(),
			   data.GetTagVNtrk(), data.GetTagVZerr(),
			   data.GetTagVChisq(), data.GetTagVNdf(),
			   data.GetTagVIsl(), dtres_param );

  //Scp component
  double sin_pdf;
      sin_pdf =
      0.5 / Tau *
      AfRkRdetRnp_fullrec( dt, 1., Tau,
			   data.GetDeltam(), data.GetAk(), data.GetCk(),
			   data.GetRecVNtrk(), data.GetRecVZerr(),
			   data.GetRecVChisq(), data.GetRecVNdf(),
			   data.GetTagVNtrk(), data.GetTagVZerr(),
			   data.GetTagVChisq(), data.GetTagVNdf(),
			   data.GetTagVIsl(), dtres_param );
	
  //pipi PDF
  double pipi_pdf = life_pdf;

  //Normalise over q
  pipi_pdf /= 2.0*int_life_pdf;

  AddOutlier(data, dt, pipi_pdf);

   return pipi_pdf;
}




void AddOutlier( const Data& data, const double& dt,
			   double& pdf )
{
  const dtres_param_t* const dtres_param =
    get_dtres_param( data.GetExpNo(), data.GetMC() );

  const double out_pdf =
    gaussian( dt, 0.0, dtres_param->sig_ol ) /
    (2.0*norm_gaussian( h_dt_ll, h_dt_ul, 0.0, dtres_param->sig_ol ));

  double f_ol;
  if( 1 < data.GetRecVNtrk() && 1 < data.GetTagVNtrk() )
    f_ol = dtres_param->fol_mul;
  else
    f_ol = dtres_param->fol_sgl;

  pdf = ( (1.0 - f_ol)*pdf ) + ( f_ol*out_pdf );

  return;
}


double QQbarDtQC( const Data& data, const double& dt,
				const int& q, 
				const std::vector<double>& par,
				unsigned int par0 )
{

//  if( data.GetExpNo() > 29 )
//    par0 += npar_qqbar_svd1_;

  const double sigma =
    sqrt( (data.GetRecVZerr()*data.GetRecVZerr()) +
	  (data.GetTagVZerr()*data.GetTagVZerr()) );

  double tau_eff, f_prompt/*, DC*/;
//  if( data.GetExpNo() < 29 )
//    {
//      tau_eff  = par[par0+0+npar_qqbar_svd1_];
//      f_prompt = par[par0+2+npar_qqbar_svd1_];
//    }
//  else
//    {
      tau_eff  = par[par0+0];
      f_prompt = par[par0+2];
//    }
  const double mean_prompt = par[par0+1];
  double mean_res, sigma_main, f_tail, sigma_tail;

  mean_res = par[par0+3];
  sigma_main = sigma*par[par0+4]*dt_resol_global::inv_bgc;
  f_tail = par[par0+5];
  sigma_tail = sigma_main*par[par0+6];

  //Lifetime component
  const double life_pdf =
    (1.0 - f_prompt) *
    ( ((1.0 - f_tail)*
       Ef_conv_gauss(dt, tau_eff, mean_res, sigma_main)) +
      (f_tail*
       Ef_conv_gauss(dt, tau_eff, mean_res, sigma_tail)) );

  const double int_life_pdf =
    (1.0 - f_prompt) *
    ( ((1.0 - f_tail)*
       norm_Ef_conv_gauss(h_dt_ll, h_dt_ul, tau_eff, mean_res, sigma_main)) +
      (f_tail*
       norm_Ef_conv_gauss(h_dt_ll, h_dt_ul, tau_eff, mean_res, sigma_tail)) );

  //Prompt component
  const double prompt_pdf =
    f_prompt *
    ( ((1.0 - f_tail)*gaussian(dt, mean_prompt, sigma_main)) +
      (f_tail*gaussian(dt, mean_prompt, sigma_tail)) );

  const double int_prompt_pdf =
    f_prompt *
    ( ((1.0 - f_tail)*
       norm_gaussian(h_dt_ll, h_dt_ul, mean_prompt, sigma_main)) +
      (f_tail*norm_gaussian(h_dt_ll, h_dt_ul, mean_prompt, sigma_tail)) );

  //QQbar PDF
  double qqbar_pdf = life_pdf + prompt_pdf;

  //Normalise over q, c
  qqbar_pdf /= 2.0*(int_life_pdf + int_prompt_pdf);

  AddOutlier(data, dt, qqbar_pdf);

  return qqbar_pdf;
}


void AddOutDtQC( const Data& data, const double& dt,
			       double& pdf )
{
  const dtres_param_t* const dtres_param =
    get_dtres_param( data.GetExpNo(), 0 );

  const double out_pdf =
    gaussian( dt, 0.0, dtres_param->sig_ol ) /
    (4.0*norm_gaussian( h_dt_ll, h_dt_ul, 0.0, dtres_param->sig_ol ));

  double f_ol;
  if( 1 < data.GetRecVNtrk() && 1 < data.GetTagVNtrk() )
    f_ol = dtres_param->fol_mul;
  else
    f_ol = dtres_param->fol_sgl;

  pdf = ( (1.0 - f_ol)*pdf ) + ( f_ol*out_pdf );

  return;
}


#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
