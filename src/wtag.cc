#include "wtag.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

//Wrong tag fraction
const int set_rbin(const double r)
{
  return(0.   <=r && r<=0.1   ? 0 :
	 0.1  < r && r<=0.25  ? 1 :
	 0.25 < r && r<=0.5   ? 2 :
	 0.5  < r && r<=0.625 ? 3 :
	 0.625< r && r<=0.75  ? 4 :
	 0.75 < r && r<=0.875 ? 5 :
	 0.875< r && r<=1.0   ? 6 : 7);
}

double set_r( const unsigned int& r_bin )
{
  if( r_bin == 0 )
    return 0.05;
  else if( r_bin == 1 )
    return 0.2;
  else if( r_bin == 2 )
    return 0.3;
  else if( r_bin == 3 )
    return 0.6;
  else if( r_bin == 4 )
    return 0.7;
  else if( r_bin == 5 )
    return 0.8;
  else
    return 0.9;
}

double set_wtag( const unsigned int& exp_no, const int& rbin,
		 const unsigned int& mc_type )
{
  const double w_mc_svd1[7] =
    { 0.5, 0.4208, 0.3003, 0.2193, 0.1546, 0.0916, 0.0229 };

  const double w_mc_svd2[7] =
    { 0.5, 0.4122, 0.3078, 0.2128, 0.1499, 0.0913, 0.0219 };

  const double w_data_svd1[7] =
    { 0.5, 0.4189, 0.3299, 0.2339, 0.1706, 0.0979, 0.0229 };

  const double w_data_svd2[7] =
    { 0.5, 0.4188, 0.3193, 0.2229, 0.1632, 0.1041, 0.0251 };

  if( exp_no < 29 && mc_type == 0 )
    return w_data_svd1[rbin];
  else if( exp_no < 29 && mc_type == 1 )
    return w_mc_svd1[rbin];
  else if( exp_no > 29 && mc_type == 0 )
    return w_data_svd2[rbin];
  else
    return w_mc_svd2[rbin];
}

double set_dwtag( const unsigned int& exp_no, const int& rbin,
		  const unsigned int& mc_type )
{
  const double dw_mc_svd1[7] =
    { 0., +0.0583, +0.0057, -0.0393, +0.0047, -0.0119, -0.0059 };

  const double dw_mc_svd2[7] =
    { 0., +0.0041, +0.0103, -0.0048, +0.0015, +0.0144, +0.0019 };

  const double dw_data_svd1[7] =
    { 0., +0.0570, +0.0126, -0.0148, -0.0006, -0.0089, -0.0047 };

  const double dw_data_svd2[7] =
    { 0., +0.0088, -0.0104, -0.0109, -0.0186, +0.0017, -0.0036 };

  if( exp_no < 29 && mc_type == 0 )
    return dw_data_svd1[rbin];
  else if( exp_no < 29 && mc_type == 1 )
    return dw_mc_svd1[rbin];
  else if( exp_no > 29 && mc_type == 0 )
    return dw_data_svd2[rbin];
  else
    return dw_mc_svd2[rbin];
}

double set_wtag_err( const unsigned int& exp_no, const int& rbin )
{
  const double w_data_svd1_err[7] =
    { 0.0, 0.0058, 0.0060, 0.0066, 0.0057, 0.0058, 0.0034 };

  const double w_data_svd2_err[7] =
    { 0.0, 0.0026, 0.0023, 0.0026, 0.0026, 0.0025, 0.0014 };

  if( exp_no < 29 )
    return w_data_svd1_err[rbin];
  else
    return w_data_svd2_err[rbin];
}

double set_dwtag_err( const unsigned int& exp_no, const int& rbin )
{
  const double dw_data_svd1_err[7] =
    { 0.0, 0.0086, 0.0090, 0.0099, 0.0088, 0.0090, 0.0056 };

  const double dw_data_svd2_err[7] =
    { 0.0, 0.0039, 0.0034, 0.0040, 0.0041, 0.0049, 0.0023 };

  if( exp_no < 29 )
    return dw_data_svd1_err[rbin];
  else
    return dw_data_svd2_err[rbin];
}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
