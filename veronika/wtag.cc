#include <iostream>
#include <fstream>
#include <sstream>

#include "tatami/tatami.h"

#include "Minuit2/MnUserParameters.h"

#include "TDirectory.h"

#include "tools.h"
#include "hist.h"
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

void load_fdhist_rbin( std::vector<Data>& data,
		       std::vector<TH1D*>& h_fd_svd1,
		       std::vector<TH1D*>& h_fd_svd2,
		       const unsigned int& id )
{
  TH1D *h_fd_rbin0_svd1 =
    new TH1D( Form("h_fd_rbin0_svd1%d",id), "", bins_fd, h_fd_ll, h_fd_ul );
  TH1D *h_fd_rbin1_svd1 =
    new TH1D( Form("h_fd_rbin1_svd1%d",id), "", bins_fd, h_fd_ll, h_fd_ul );
  TH1D *h_fd_rbin2_svd1 =
    new TH1D( Form("h_fd_rbin2_svd1%d",id), "", bins_fd, h_fd_ll, h_fd_ul );
  TH1D *h_fd_rbin3_svd1 =
    new TH1D( Form("h_fd_rbin3_svd1%d",id), "", bins_fd, h_fd_ll, h_fd_ul );
  TH1D *h_fd_rbin4_svd1 =
    new TH1D( Form("h_fd_rbin4_svd1%d",id), "", bins_fd, h_fd_ll, h_fd_ul );
  TH1D *h_fd_rbin5_svd1 =
    new TH1D( Form("h_fd_rbin5_svd1%d",id), "", bins_fd, h_fd_ll, h_fd_ul );
  TH1D *h_fd_rbin6_svd1 =
    new TH1D( Form("h_fd_rbin6_svd1%d",id), "", bins_fd, h_fd_ll, h_fd_ul );

  TH1D *h_fd_rbin0_svd2 =
    new TH1D( Form("h_fd_rbin0_svd2%d",id), "", bins_fd, h_fd_ll, h_fd_ul );
  TH1D *h_fd_rbin1_svd2 =
    new TH1D( Form("h_fd_rbin1_svd2%d",id), "", bins_fd, h_fd_ll, h_fd_ul );
  TH1D *h_fd_rbin2_svd2 =
    new TH1D( Form("h_fd_rbin2_svd2%d",id), "", bins_fd, h_fd_ll, h_fd_ul );
  TH1D *h_fd_rbin3_svd2 =
    new TH1D( Form("h_fd_rbin3_svd2%d",id), "", bins_fd, h_fd_ll, h_fd_ul );
  TH1D *h_fd_rbin4_svd2 =
    new TH1D( Form("h_fd_rbin4_svd2%d",id), "", bins_fd, h_fd_ll, h_fd_ul );
  TH1D *h_fd_rbin5_svd2 =
    new TH1D( Form("h_fd_rbin5_svd2%d",id), "", bins_fd, h_fd_ll, h_fd_ul );
  TH1D *h_fd_rbin6_svd2 =
    new TH1D( Form("h_fd_rbin6_svd2%d",id), "", bins_fd, h_fd_ll, h_fd_ul );

  for( std::vector<Data>::iterator it = data.begin();
       it != data.end(); it++ )
    {
      Data data_ = *(it);

      const int r_bin = set_rbin( data_.GetR() );

      if( r_bin == 0 )
	if( data_.GetExpNo() < 29 )
	  h_fd_rbin0_svd1->Fill( data_.GetFD() );
	else
	  h_fd_rbin0_svd2->Fill( data_.GetFD() );
      else if( r_bin == 1 )
	if( data_.GetExpNo() < 29 )
	  h_fd_rbin1_svd1->Fill( data_.GetFD() );
	else
	  h_fd_rbin1_svd2->Fill( data_.GetFD() );
      else if( r_bin == 2 )
	if( data_.GetExpNo() < 29 )
	  h_fd_rbin2_svd1->Fill( data_.GetFD() );
	else
	  h_fd_rbin2_svd2->Fill( data_.GetFD() );
      else if( r_bin == 3 )
	if( data_.GetExpNo() < 29 )
	  h_fd_rbin3_svd1->Fill( data_.GetFD() );
	else
	  h_fd_rbin3_svd2->Fill( data_.GetFD() );
      else if( r_bin == 4 )
	if( data_.GetExpNo() < 29 )
	  h_fd_rbin4_svd1->Fill( data_.GetFD() );
	else
	  h_fd_rbin4_svd2->Fill( data_.GetFD() );
      else if( r_bin == 5 )
	if( data_.GetExpNo() < 29 )
	  h_fd_rbin5_svd1->Fill( data_.GetFD() );
	else
	  h_fd_rbin5_svd2->Fill( data_.GetFD() );
      else
	if( data_.GetExpNo() < 29 )
	  h_fd_rbin6_svd1->Fill( data_.GetFD() );
	else
	  h_fd_rbin6_svd2->Fill( data_.GetFD() );
    }
  h_fd_rbin0_svd1->Sumw2();
  h_fd_rbin1_svd1->Sumw2();
  h_fd_rbin2_svd1->Sumw2();
  h_fd_rbin3_svd1->Sumw2();
  h_fd_rbin4_svd1->Sumw2();
  h_fd_rbin5_svd1->Sumw2();
  h_fd_rbin6_svd1->Sumw2();
  h_fd_rbin0_svd2->Sumw2();
  h_fd_rbin1_svd2->Sumw2();
  h_fd_rbin2_svd2->Sumw2();
  h_fd_rbin3_svd2->Sumw2();
  h_fd_rbin4_svd2->Sumw2();
  h_fd_rbin5_svd2->Sumw2();
  h_fd_rbin6_svd2->Sumw2();
  h_fd_rbin0_svd1->Scale(1.0/norm_hist1d(h_fd_ll, h_fd_ul,
					 h_fd_rbin0_svd1));
  h_fd_rbin1_svd1->Scale(1.0/norm_hist1d(h_fd_ll, h_fd_ul,
					 h_fd_rbin1_svd1));
  h_fd_rbin2_svd1->Scale(1.0/norm_hist1d(h_fd_ll, h_fd_ul,
					 h_fd_rbin2_svd1));
  h_fd_rbin3_svd1->Scale(1.0/norm_hist1d(h_fd_ll, h_fd_ul,
					 h_fd_rbin3_svd1));
  h_fd_rbin4_svd1->Scale(1.0/norm_hist1d(h_fd_ll, h_fd_ul,
					 h_fd_rbin4_svd1));
  h_fd_rbin5_svd1->Scale(1.0/norm_hist1d(h_fd_ll, h_fd_ul,
					 h_fd_rbin5_svd1));
  h_fd_rbin6_svd1->Scale(1.0/norm_hist1d(h_fd_ll, h_fd_ul,
					 h_fd_rbin6_svd1));
  h_fd_rbin0_svd2->Scale(1.0/norm_hist1d(h_fd_ll, h_fd_ul,
					 h_fd_rbin0_svd2));
  h_fd_rbin1_svd2->Scale(1.0/norm_hist1d(h_fd_ll, h_fd_ul,
					 h_fd_rbin1_svd2));
  h_fd_rbin2_svd2->Scale(1.0/norm_hist1d(h_fd_ll, h_fd_ul,
					 h_fd_rbin2_svd2));
  h_fd_rbin3_svd2->Scale(1.0/norm_hist1d(h_fd_ll, h_fd_ul,
					 h_fd_rbin3_svd2));
  h_fd_rbin4_svd2->Scale(1.0/norm_hist1d(h_fd_ll, h_fd_ul,
					 h_fd_rbin4_svd2));
  h_fd_rbin5_svd2->Scale(1.0/norm_hist1d(h_fd_ll, h_fd_ul,
					 h_fd_rbin5_svd2));
  h_fd_rbin6_svd2->Scale(1.0/norm_hist1d(h_fd_ll, h_fd_ul,
					 h_fd_rbin6_svd2));

  h_fd_svd1.push_back(h_fd_rbin0_svd1);
  h_fd_svd1.push_back(h_fd_rbin1_svd1);
  h_fd_svd1.push_back(h_fd_rbin2_svd1);
  h_fd_svd1.push_back(h_fd_rbin3_svd1);
  h_fd_svd1.push_back(h_fd_rbin4_svd1);
  h_fd_svd1.push_back(h_fd_rbin5_svd1);
  h_fd_svd1.push_back(h_fd_rbin6_svd1);
  h_fd_svd2.push_back(h_fd_rbin0_svd2);
  h_fd_svd2.push_back(h_fd_rbin1_svd2);
  h_fd_svd2.push_back(h_fd_rbin2_svd2);
  h_fd_svd2.push_back(h_fd_rbin3_svd2);
  h_fd_svd2.push_back(h_fd_rbin4_svd2);
  h_fd_svd2.push_back(h_fd_rbin5_svd2);
  h_fd_svd2.push_back(h_fd_rbin6_svd2);
}

double set_wtag( const unsigned int& exp_no, const unsigned int& mc_type,
		 const int& rbin )
{
  const double w_svd1_mc[7] = { 0.5, 0.413653, 0.310305, 0.215106, 0.151385,
				0.0934167, 0.0212698 };

  const double w_svd2_mc[7] = { 0.5, 0.412019, 0.304503, 0.208394, 0.149832,
				0.0865786, 0.0243203 };

  const double w_svd1_data[7] = { 0.5, 0.422837, 0.336574, 0.235379, 0.166249,
				  0.104922, 0.0262446 };

  const double w_svd2_data[7] = { 0.5, 0.429423, 0.327273, 0.22304, 0.160854,
				  0.105316, 0.0192838 };

  if( exp_no < 29 && mc_type == 0 )
    return w_svd1_data[rbin];
  else if( exp_no < 29 && mc_type == 1 )
    return w_svd1_mc[rbin];
  else if( exp_no > 29 && mc_type == 0 )
    return w_svd2_data[rbin];
  else
    return w_svd2_mc[rbin];
}

double set_dwtag( const unsigned int& exp_no, const unsigned int& mc_type,
		  const int& rbin )
{
  const double dw_svd1_mc[7] = { 0., 0.0456694, -0.00811093, -0.0234073,
				 -0.00217406, -0.00573868, -0.00115232};

  const double dw_svd2_mc[7] = { 0., -0.0322979, -0.0255096, 0.0190514,
				 0.00594418, -0.0171155, 0.00400088};

  const double dw_svd1_data[7] = { 0., 0.0577827, 0.0124391, -0.0122602,
				   -0.0108194, 0.00818429, 0.0034784 };

  const double dw_svd2_data[7] = { 0., -0.039313, -0.035766, 0.0175526,
				   0.00232484, -0.0270618, -0.000737663 };

  if( exp_no < 29 && mc_type == 0 )
    return dw_svd1_data[rbin];
  else if( exp_no < 29 && mc_type == 1 )
    return dw_svd1_mc[rbin];
  else if( exp_no > 29 && mc_type == 0 )
    return dw_svd2_data[rbin];
  else
    return dw_svd2_mc[rbin];
}

const double set_wtag_uerr(const int expno, const int rbin)
{
  static const double w_svd1_uerr[7]=
    {0.0, 7.522044e-03, 7.882328e-03, 1.020837e-02, 8.434339e-03,
     7.855777e-03, 5.610706e-03};

  static const double w_svd2_uerr[7]=
    {0.0, 5.033855e-03, 6.072289e-03, 1.078750e-02, 1.039213e-02,
     6.848312e-03, 4.011759e-03};

  if(expno < 30){
    return w_svd1_uerr[rbin];
  }else{
    return w_svd2_uerr[rbin];
  }
}

const double set_wtag_lerr(const int expno, const int rbin)
{
  static const double w_svd1_lerr[7]=
    {0.0, 6.665670e-03, 8.470807e-03, 7.976266e-03, 6.813315e-03,
     6.870697e-03, 5.401010e-03};

  static const double w_svd2_lerr[7]=
    {0.0, 4.713851e-03, 5.700081e-03, 5.819431e-03, 6.491131e-03,
     7.627830e-03, 5.163255e-03};

  if(expno < 30){
    return w_svd1_lerr[rbin];
  }else{
    return w_svd2_lerr[rbin];
  }
}

const double set_delta_wtag_uerr(const int expno, const int rbin)
{
  static const double dw_svd1_uerr[7]=
    {0.0, 9.545226e-03, 9.795511e-03, 1.029756e-02, 9.826752e-03,
     9.381480e-03, 5.739762e-03};

  static const double dw_svd2_uerr[7]=
    {0.0, 6.061285e-03, 6.141136e-03, 6.605191e-03, 6.007190e-03,
     6.183453e-03, 3.652403e-03};

  if(expno < 30){
    return dw_svd1_uerr[rbin];
  }else{
    return dw_svd2_uerr[rbin];
  }
}

const double set_delta_wtag_lerr(const int expno, const int rbin)
{
  static const double dw_svd1_lerr[7]=
    {0.0, 9.454471e-03, 9.894429e-03, 1.027626e-02, 8.881274e-03,
     9.159412e-03, 5.761282e-03};

  static const double dw_svd2_lerr[7]=
    {0.0, 6.619413e-03, 6.041260e-03, 7.192499e-03, 6.075392e-03,
     5.925455e-03, 3.881354e-03};

  if(expno < 30){
    return dw_svd1_lerr[rbin];
  }else{
    return dw_svd2_lerr[rbin];
  }
}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
