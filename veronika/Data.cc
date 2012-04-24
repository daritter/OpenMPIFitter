#include<cstdlib>
#include<iostream>
#include<fstream>
#include<cmath>

#include "tatami/tatami.h"

#include "TVectorD.h"

#include "Data.h"
#include "Fisher_svd1.C"
#include "Fisher_svd2.C"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

void FillPlot( std::ifstream& datafile1, std::ifstream& datafile2, std::vector<Data>& data,
	       const unsigned int& mc_type, const unsigned int& B_type,
	       const unsigned int& mc_keep,
	       const double& ctbto_ul,
	       const double& mbc_ll, const double& mbc_ul,
	       const double& de_ll,  const double& de_ul,
	       const double& homega_ll, const double& homega_ul )
{
  /* mc_keep
     0: All
     1: Correctly reconstructed
     2: Mis-reconstructed */

  if( mc_keep != 0 && mc_keep != 1 && mc_keep != 2 )
    {
      std::cout << "ERROR: invalid mc_keep. Check code." << std::endl;
      exit(1);
    }

  if( datafile1.is_open() == 0 )
    {
      std::cerr << "ERROR: Datafile could not be opened."
                << std::endl;
      exit(1);
    }

  //double dummy; //Represents input value

  double dummy1_exp_no;
  double dummy2_run_no;
  double dummy3_evt_no;

  double dummy4_mbc;
  double dummy5_benergy;
  double dummy6_de;
  double dummy7_lr;
  double dummy8_momega;
  double dummy9_homega;
  double dummy10_background;

  double dummy11_good_mc;
  double dummy12_momega_gen;

  double dummy13_c;

  double dummy14_delta_z;
  double dummy15_cpv_chisq;
  double dummy16_cpv_ndf;
  double dummy17_cpv_zerr;
  double dummy18_cpv_ntrk;

  double dummy19_q;
  double dummy20_r;

  double dummy21_tgv_chisq;
  double dummy22_tgv_ndf;
  double dummy23_tgv_zerr;
  double dummy24_tgv_ntrk;
  double dummy25_tgv_isl;

  double dummy26_costheta_b;

  double dummy40_l0c;
  double dummy41_l2c;
  double dummy42_l0n;
  double dummy43_l2n;
  double dummy44_ctbto;
  double dummy45_ctbz;
  double dummy46_cbz;
  double dummy47_sumpt;

  datafile1 >> dummy1_exp_no >> dummy2_run_no >> dummy3_evt_no
	    >> dummy4_mbc >> dummy5_benergy >> dummy6_de
	    >> dummy7_lr >> dummy8_momega >> dummy9_homega >> dummy10_background
	    >> dummy11_good_mc >> dummy12_momega_gen
	    >> dummy13_c >> dummy14_delta_z
	    >> dummy15_cpv_chisq >> dummy16_cpv_ndf
	    >> dummy17_cpv_zerr >> dummy18_cpv_ntrk
	    >> dummy19_q >> dummy20_r
	    >> dummy21_tgv_chisq >> dummy22_tgv_ndf >> dummy23_tgv_zerr
	    >> dummy24_tgv_ntrk >> dummy25_tgv_isl >> dummy26_costheta_b;

  datafile2 >> dummy40_l0c >> dummy41_l2c
	    >> dummy42_l0n >> dummy43_l2n
	    >> dummy44_ctbto >> dummy45_ctbz
	    >> dummy46_cbz >> dummy47_sumpt;
  while( datafile1.eof() == 0 )
    {
      Data data_;

      if( //dummy44_ctbto <= ctbto_ul &&
	  mbc_ll   <= dummy4_mbc /*&& dummy4_mbc <= mbc_ul*/   &&
      	  de_ll    <= dummy6_de  && dummy6_de  <= de_ul    &&
      	  h_lr_ll <= dummy7_lr &&
	  h_momega_ll <= dummy8_momega && dummy8_momega <= h_momega_ul &&
	  //h_homega_ll   <= dummy9_homega && dummy9_homega <= h_homega_ul   &&
	  fabs(dummy14_delta_z*Belle::dt_resol_global::inv_bgc) < dt_cut &&
	  (dummy16_cpv_ndf == 0.0 ||
	   dummy15_cpv_chisq/dummy16_cpv_ndf < h_cut) &&
	  (dummy22_tgv_ndf == 0.0 ||
	   dummy21_tgv_chisq/dummy22_tgv_ndf < h_cut) &&
	  ((dummy16_cpv_ndf == 0.0 && dummy17_cpv_zerr < sigmaz_sngl_cut) ||
	   (dummy16_cpv_ndf >  0.0 && dummy17_cpv_zerr < sigmaz_mult_cut)) &&
	  ((dummy22_tgv_ndf == 0.0 && dummy23_tgv_zerr < sigmaz_sngl_cut) ||
	   (dummy22_tgv_ndf >  0.0 && dummy23_tgv_zerr < sigmaz_mult_cut)) 
	   )
	{
	  data_.SetExpNo( static_cast<unsigned int>(dummy1_exp_no) );
	  data_.SetRunNo( static_cast<unsigned int>(dummy2_run_no) );
	  data_.SetEvtNo( static_cast<unsigned int>(dummy3_evt_no) );

	  data_.SetMbc( dummy4_mbc );
	  data_.SetBenergy( dummy5_benergy );
	  data_.SetDE( dummy6_de );
	  data_.SetLR( dummy7_lr);
	  data_.SetFD( log( (dummy7_lr-h_lr_ll)/(h_lr_ul-dummy7_lr) ) );

	  data_.SetMomega( dummy8_momega );
	  data_.SetHomega( dummy9_homega );

	  if(static_cast<unsigned int>(dummy10_background) == 1)
	  data_.SetGoodMC(true);
	  else
	  data_.SetGoodMC(false);
	  data_.SetMomegaGen( dummy12_momega_gen );

	  data_.SetMC(mc_type);
	  data_.SetBtype(B_type);
          if( data_.GetMC() == 1 && data_.GetBtype() == 0 )
	    {
	      data_.SetTau(b0_lifetime_mc);
              data_.SetDeltam(delta_m_mc);
	    }
          if( data_.GetMC() == 0 && data_.GetBtype() == 0 )
	    {
	      data_.SetTau(b0_lifetime);
              data_.SetDeltam(delta_m);
	    }

	  data_.SetDeltat(dummy14_delta_z*Belle::dt_resol_global::inv_bgc);

	  data_.SetRecVChisq( dummy15_cpv_chisq );
	  data_.SetRecVNdf( static_cast<unsigned int>(dummy16_cpv_ndf) );
	  data_.SetRecVZerr( dummy17_cpv_zerr );
	  data_.SetRecVNtrk( static_cast<unsigned int>(dummy18_cpv_ntrk) );

          data_.SetQ( static_cast<int>(dummy19_q) );
          data_.SetR( fabs(dummy20_r) );

	  if( dummy13_c == 1.0 )
	    data_.SetC(1);
	  if( dummy13_c == 2.0 )
	    data_.SetC(-1);

	  data_.SetTagVChisq( dummy21_tgv_chisq );
	  data_.SetTagVNdf( static_cast<unsigned int>(dummy22_tgv_ndf) );
	  data_.SetTagVZerr( dummy23_tgv_zerr );
	  data_.SetTagVNtrk( static_cast<unsigned int>(dummy24_tgv_ntrk) );
	  data_.SetTagVIsl(static_cast<int>(dummy25_tgv_isl));

	  double dummy_a_k = 0.0;
	  double dummy_c_k = 0.0;
	  const double pb_cms_sq =
	    (dummy5_benergy*dummy5_benergy) - (dummy4_mbc*dummy4_mbc);
	  const double Eb_cms =
	    sqrt( (Belle::dt_resol_global::mbzero*
		   Belle::dt_resol_global::mbzero) +
		  pb_cms_sq );
	  Belle::CalcAkCk( dummy26_costheta_b, Eb_cms, &dummy_a_k, &dummy_c_k,
			   Belle::dt_resol_global::mbzero );
	  data_.SetAk( dummy_a_k );
	  data_.SetCk( dummy_c_k );
	  data_.SetCosTheta( dummy26_costheta_b );

	  //Resolution function calculation
	  const dtres_param_t* const dtres_param =
	    get_dtres_param( data_.GetExpNo(), data_.GetMC() );
	  const double life_pdf =
	    EfRkRdetRnp_fullrec( data_.GetDeltat(), 
				 data_.GetBtype(),  data_.GetTau(),
				 data_.GetAk(),        data_.GetCk(),
				 data_.GetRecVNtrk(),  data_.GetRecVZerr(),
				 data_.GetRecVChisq(), data_.GetRecVNdf(),
				 data_.GetTagVNtrk(),  data_.GetTagVZerr(),
				 data_.GetTagVChisq(), data_.GetTagVNdf(),
				 data_.GetTagVIsl(),   dtres_param );
	  data_.SetEfRkRdetRnp( life_pdf );
	  const double int_life_pdf =
	    norm_EfRkRdetRnp_fullrec( dt_resol_global::dt_llmt,
				      dt_resol_global::dt_ulmt,
				      data_.GetBtype(),  data_.GetTau(),
				      data_.GetAk(),        data_.GetCk(),
				      data_.GetRecVNtrk(),  data_.GetRecVZerr(),
				      data_.GetRecVChisq(), data_.GetRecVNdf(),
				      data_.GetTagVNtrk(),  data_.GetTagVZerr(),
				      data_.GetTagVChisq(), data_.GetTagVNdf(),
				      data_.GetTagVIsl(),   dtres_param );
	  data_.SetNormEfRkRdetRnp( int_life_pdf );
	  const double cos_pdf =
	    0.5 / data_.GetTau() *
	    MfRkRdetRnp_fullrec( data_.GetDeltat(), 
				 data_.GetBtype(),  data_.GetTau(),
				 data_.GetDeltam(),
				 data_.GetAk(),        data_.GetCk(),
				 data_.GetRecVNtrk(),  data_.GetRecVZerr(),
				 data_.GetRecVChisq(), data_.GetRecVNdf(),
				 data_.GetTagVNtrk(),  data_.GetTagVZerr(),
				 data_.GetTagVChisq(), data_.GetTagVNdf(),
				 data_.GetTagVIsl(),   dtres_param );
	  data_.SetMfRkRdetRnp( cos_pdf );
	  const double sin_pdf =
	    0.5 / data_.GetTau() *
	    AfRkRdetRnp_fullrec( data_.GetDeltat(), 
				 data_.GetBtype(),  data_.GetTau(),
				 data_.GetDeltam(),
				 data_.GetAk(),        data_.GetCk(),
				 data_.GetRecVNtrk(),  data_.GetRecVZerr(),
				 data_.GetRecVChisq(), data_.GetRecVNdf(),
				 data_.GetTagVNtrk(),  data_.GetTagVZerr(),
				 data_.GetTagVChisq(), data_.GetTagVNdf(),
				 data_.GetTagVIsl(),   dtres_param );
	  data_.SetAfRkRdetRnp( sin_pdf );

	  if( h_fd_ll < data_.GetFD() && data_.GetFD() < h_fd_ul )
	    switch( mc_keep )
	      {
	      case 1:
		if( data_.GetGoodMC() == true )
		  data.push_back(data_);
		break;

	      case 2:
		if( data_.GetGoodMC() == false )
		  data.push_back(data_);
		break;

	      default:
		data.push_back(data_);
		break;
	      }
	}

      datafile1 >> dummy1_exp_no >> dummy2_run_no >> dummy3_evt_no
		>> dummy4_mbc >> dummy5_benergy >> dummy6_de
		>> dummy7_lr >> dummy8_momega >> dummy9_homega >> dummy10_background
	        >> dummy11_good_mc >> dummy12_momega_gen
		>> dummy13_c >> dummy14_delta_z
		>> dummy15_cpv_chisq >> dummy16_cpv_ndf
		>> dummy17_cpv_zerr >> dummy18_cpv_ntrk
		>> dummy19_q >> dummy20_r
		>> dummy21_tgv_chisq >> dummy22_tgv_ndf >> dummy23_tgv_zerr
		>> dummy24_tgv_ntrk >> dummy25_tgv_isl >> dummy26_costheta_b;

      datafile2 >> dummy40_l0c >> dummy41_l2c
		>> dummy42_l0n >> dummy43_l2n
		>> dummy44_ctbto >> dummy45_ctbz
		>> dummy46_cbz >> dummy47_sumpt;
    }
  datafile1.close();
  datafile2.close();
  std::cout << data.size() << " events read." 
	    << std::endl << std::endl;
}

void FillPlotGen( std::ifstream& datafile, std::vector<Data>& data )
{
  if( datafile.is_open() == 0 )
    {
      std::cerr << "ERROR: Datafile could not be opened."
                << std::endl;
      exit(1);
    }

  double dummy12_momega_gen; //Represents input value

  datafile >> dummy12_momega_gen;
  while( datafile.eof() == 0 )
    {
      Data data_;

      if( h_momega_ll < dummy12_momega_gen && dummy12_momega_gen < h_momega_ul )
	{
	  data_.SetMomega( dummy12_momega_gen );
	  data_.SetMomegaGen( dummy12_momega_gen );

	  data.push_back(data_);
	}

      datafile >> dummy12_momega_gen;
    }
  datafile.close();
  std::cout << data.size() << " events read." 
	    << std::endl << std::endl;
}

std::vector<unsigned int> strip_expno( std::vector<Data>& data )
{
  std::vector<unsigned int> strip;
  for( std::vector<Data>::iterator it = data.begin();
       it != data.end(); it++ )
    {
      Data data_ = *(it);

      strip.push_back( data_.GetExpNo() );
    }
  return strip;
}

std::vector<double> strip_recv_chisq( std::vector<Data>& data )
{
  std::vector<double> strip;
  for( std::vector<Data>::iterator it = data.begin();
       it != data.end(); it++ )
    {
      Data data_ = *(it);

      strip.push_back( data_.GetRecVChisq() );
    }
  return strip;
}

std::vector<unsigned int> strip_recv_ndf( std::vector<Data>& data )
{
  std::vector<unsigned int> strip;
  for( std::vector<Data>::iterator it = data.begin();
       it != data.end(); it++ )
    {
      Data data_ = *(it);

      strip.push_back( data_.GetRecVNdf() );
    }
  return strip;
}

std::vector<double> strip_recv_zerr( std::vector<Data>& data )
{
  std::vector<double> strip;
  for( std::vector<Data>::iterator it = data.begin();
       it != data.end(); it++ )
    {
      Data data_ = *(it);

      strip.push_back( data_.GetRecVZerr() );
    }
  return strip;
}

std::vector<double> strip_tagv_chisq( std::vector<Data>& data )
{
  std::vector<double> strip;
  for( std::vector<Data>::iterator it = data.begin();
       it != data.end(); it++ )
    {
      Data data_ = *(it);

      strip.push_back( data_.GetTagVChisq() );
    }
  return strip;
}

std::vector<unsigned int> strip_tagv_ndf( std::vector<Data>& data )
{
  std::vector<unsigned int> strip;
  for( std::vector<Data>::iterator it = data.begin();
       it != data.end(); it++ )
    {
      Data data_ = *(it);

      strip.push_back( data_.GetTagVNdf() );
    }
  return strip;
}

std::vector<double> strip_tagv_zerr( std::vector<Data>& data )
{
  std::vector<double> strip;
  for( std::vector<Data>::iterator it = data.begin();
       it != data.end(); it++ )
    {
      Data data_ = *(it);

      strip.push_back( data_.GetTagVZerr() );
    }
  return strip;
}

std::vector<unsigned int> strip_tagv_isl( std::vector<Data>& data )
{
  std::vector<unsigned int> strip;
  for( std::vector<Data>::iterator it = data.begin();
       it != data.end(); it++ )
    {
      Data data_ = *(it);

      strip.push_back( data_.GetTagVIsl() );
    }
  return strip;
}

std::vector<Data> Add( std::vector<Data>& data1, std::vector<Data>& data2 )
{
  std::vector<Data> data;

  for( std::vector<Data>::iterator it = data1.begin();
       it != data1.end(); it++ )
    {
      Data data_ = *(it);

      data.push_back(data_);
    }

  for( std::vector<Data>::iterator it = data2.begin();
       it != data2.end(); it++ )
    {
      Data data_ = *(it);

      data.push_back(data_);
    }

  return data;
}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
