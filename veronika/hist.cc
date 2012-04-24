#include <iostream>
#include <fstream>
#include <sstream>

#include "tatami/tatami.h"

#include "Minuit2/MnUserParameters.h"

#include "TROOT.h"
#include "TDirectory.h"


#include "constant.h"
#include "hist.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

// Canvas style
void set_style()
{
  TStyle *style = new TStyle("style","");
  style->SetCanvasBorderMode(0);
  style->SetPadBorderMode(0);
  style->SetFrameBorderMode(0);
  style->SetPadColor(0);
  style->SetCanvasColor(0);
  style->SetStatColor(0);
  style->SetPalette(1,0);

  style->SetPadLeftMargin(0.25);
  style->SetPadBottomMargin(0.20);
  //style->SetOptStat(kFALSE);
  style->SetHistMinimumZero(kTRUE);
  style->SetHistLineWidth(2);
  style->SetOptStat(0);  

  style->SetLabelSize(0.07, "XY");
  style->SetTitleSize(0.08, "XY");

  style->SetLabelOffset(0.02, "X");
  style->SetTitleOffset(1.0, "X");

  style->SetTitleOffset(1.6, "Y");
  
  style->SetPadTickX(1);
  style->SetPadTickY(1);

  style->cd();
}

// Histogram style
void hist_style( TH1* hist, const std::string& title_x,
		 const std::string& units )
{
  std::string s_xtitle;
  if( units.length() == 0 )
    s_xtitle = title_x;
  else
    s_xtitle = title_x + " [" + units + "]";

  const double bin_width =
    (hist->GetXaxis()->GetXmax() - hist->GetXaxis()->GetXmin())/
    static_cast<double>(hist->GetXaxis()->GetNbins());
  std::stringstream ss_ytitle;
  ss_ytitle << bin_width;

  std::string s_ytitle;
  if( units.length() == 0 )
    s_ytitle = "Events / [" + ss_ytitle.str() + "]";
  else
    s_ytitle = "Events / [" + ss_ytitle.str() + " " + units + "]";

  hist->GetXaxis()->SetNdivisions(50205, kFALSE);
  hist->GetXaxis()->SetTitle(s_xtitle.c_str());
  hist->GetYaxis()->SetTitle(s_ytitle.c_str());
}

// 1D Histogram PDF
double hist1d( const double& x,
	       const TH1D* h_ref )
{
  const unsigned int bin_x = h_ref->GetXaxis()->FindBin(x);

  const double hist_pdf = h_ref->GetBinContent(bin_x);

  return hist_pdf;
}

double norm_hist1d( const double& x_ll, const double& x_ul,
		    const TH1D* h_ref )
{
  //Boundary effects
  const unsigned int bin_x_ll = h_ref->GetXaxis()->FindBin(x_ll+1.0e-10);
  const unsigned int bin_x_ul = h_ref->GetXaxis()->FindBin(x_ul-1.0e-10);
  
  const double norm_hist_pdf = h_ref->Integral(bin_x_ll, bin_x_ul, "width");

  return norm_hist_pdf;
}

void load_hist1d( std::ifstream& datafile, TH1D* h_ )
{
  if( datafile.is_open() == 0 )
    {
      std::cerr << "ERROR: Datafile could not be opened."
                << std::endl;
      exit(1);
    }

  std::vector<double> entries;

  double dummy;

  datafile >> dummy;
  while( datafile.eof() == 0 )
    {
      entries.push_back(dummy);
      datafile >> dummy;
    }
  datafile.close();

  if( static_cast<int>(entries.size()) != h_->GetXaxis()->GetNbins() )
    {
      std::cerr << "ERROR: Bin inconsistency."
                << std::endl;
      exit(1);
    }

  for( unsigned int i=0; i<entries.size(); i++ )
    h_->SetBinContent(i+1, entries.at(i) );
}

//Symmetrise histogram
void symm_hist1d( TH1D* h_ )
{
  for( int i=0; i<h_->GetXaxis()->GetNbins(); i++ )
    {
      const double entries =
	(h_->GetBinContent(i+1) +
	 h_->GetBinContent(h_->GetXaxis()->GetNbins()-i))/
	2.0;

      h_->SetBinContent(i+1, entries);
      h_->SetBinContent(h_->GetXaxis()->GetNbins()-i, entries);
    }
}

// Calculate chi-square
double get_chisq( const TH1* h_data, const TH1* h_pdf,
		  const unsigned int& free_par )
{
  if( h_data->GetXaxis()->GetNbins() != h_pdf->GetXaxis()->GetNbins() )
    {
      std::cout << "Error: Data and PDF histograms must have equal bins"
		<< std::endl;
      exit(1);
    }

  double chisq = 0.0;
  unsigned int bins = 0;
  for( int i=0; i<h_data->GetXaxis()->GetNbins(); i++ )
    {
      const double x_data = h_data->GetBinContent(i+1);
      const double x_pdf  = h_pdf->GetBinContent(i+1);
      const double x_err  = sqrt(x_data);

      if( x_data != 0.0 )
	{
	  chisq += (x_data-x_pdf)*(x_data-x_pdf)/x_err/x_err;
	  bins++;
	}
    }

  return chisq/static_cast<double>(bins-free_par);
}

// 2D Histogram PDF
double hist2d( const double& x, const double& y,
	       const TH2D* h_ref )
{
  const unsigned int bin_x = h_ref->GetXaxis()->FindBin(x);
  const unsigned int bin_y = h_ref->GetYaxis()->FindBin(y);

  const double hist_pdf = h_ref->GetBinContent(bin_x, bin_y);

  return hist_pdf;
}

double norm_hist2d( const double& x_ll, const double& x_ul,
		    const double& y_ll, const double& y_ul,
		    const TH2D* h_ref )
{
  //Boundary effects
  const unsigned int bin_x_ll = h_ref->GetXaxis()->FindBin(x_ll+1.0e-10);
  const unsigned int bin_x_ul = h_ref->GetXaxis()->FindBin(x_ul-1.0e-10);
  const unsigned int bin_y_ll = h_ref->GetYaxis()->FindBin(y_ll+1.0e-10);
  const unsigned int bin_y_ul = h_ref->GetYaxis()->FindBin(y_ul-1.0e-10);

  const double norm_hist_pdf =
    h_ref->Integral(bin_x_ll, bin_x_ul, bin_y_ll, bin_y_ul, "width");

  return norm_hist_pdf;
}

// 3D Histogram PDF
double hist3d( const double& x, const double& y, const double& z,
	       const TH3D* h_ref )
{
  //Evaluate histogram PDF
  double hist_pdf = 0.0;

  const unsigned int bin_x = h_ref->GetXaxis()->FindBin(x);
  const unsigned int bin_y = h_ref->GetYaxis()->FindBin(y);
  const unsigned int bin_z = h_ref->GetZaxis()->FindBin(z);

  hist_pdf = h_ref->GetBinContent(bin_x, bin_y, bin_z);

  return hist_pdf;
}

double norm_hist3d( const double& x_ll, const double& x_ul,
		    const double& y_ll, const double& y_ul,
		    const double& z_ll, const double& z_ul,
		    const TH3D* h_ref )
{
  const unsigned int bin_x_ll = h_ref->GetXaxis()->FindBin(x_ll+1.0e-10);
  const unsigned int bin_x_ul = h_ref->GetXaxis()->FindBin(x_ul-1.0e-10);
  const unsigned int bin_y_ll = h_ref->GetYaxis()->FindBin(y_ll+1.0e-10);
  const unsigned int bin_y_ul = h_ref->GetYaxis()->FindBin(y_ul-1.0e-10);
  const unsigned int bin_z_ll = h_ref->GetZaxis()->FindBin(z_ll+1.0e-10);
  const unsigned int bin_z_ul = h_ref->GetZaxis()->FindBin(z_ul-1.0e-10);

  const double norm_hist_pdf =
    h_ref->Integral(bin_x_ll, bin_x_ul, bin_y_ll, bin_y_ul, bin_z_ll, bin_z_ul,
		    "width");

  return norm_hist_pdf;
}

// Load histograms and normalise
void load_dehist( std::vector<Data>& data, TH1D* h_de_svd1, TH1D* h_de_svd2 )
{
  for( std::vector<Data>::iterator it = data.begin();
       it != data.end(); it++ )
    {
      Data data_ = *(it);

      if( data_.GetExpNo() < 29 )
	h_de_svd1->Fill( data_.GetDE() );
      else
	h_de_svd2->Fill( data_.GetDE() );
    }
  h_de_svd1->Sumw2();
  h_de_svd2->Sumw2();
  h_de_svd1->Scale(1.0/norm_hist1d(h_de_ll, h_de_ul, h_de_svd1));
  h_de_svd2->Scale(1.0/norm_hist1d(h_de_ll, h_de_ul, h_de_svd2));
}


void load_rbinhist( std::vector<Data>& data,
		    TH1D* h_rbin_svd1, TH1D* h_rbin_svd2 )
{
  for( std::vector<Data>::iterator it = data.begin();
       it != data.end(); it++ )
    {
      Data data_ = *(it);

      double dummy_r = data_.GetR();
      if( dummy_r == 1.0 )
	dummy_r = 0.99; //For histogram boundary effect

      if( data_.GetExpNo() < 29 )
	h_rbin_svd1->Fill( dummy_r );
      else
	h_rbin_svd2->Fill( dummy_r );
    }
}

void normres_style( TCanvas* canvas )
{
  canvas->Divide(2,1);
  canvas->cd(1)->SetPad(0,0.33,1,1);
  canvas->cd(2)->SetPad(0,0,1,0.33);

  canvas->cd(1)->SetBottomMargin(0);
  canvas->cd(2)->SetTopMargin(0);
  canvas->cd(2)->SetBottomMargin(0.35);
}

void normres_style( TH1* hist, std::vector<TLine>& normres_sigma,
		    const std::string& title_x, const std::string& units )
{
  hist_style(hist, title_x, units);

  hist->GetXaxis()->SetLabelSize(0.15);
  hist->GetXaxis()->SetTitleOffset(1.0);
  hist->GetXaxis()->SetTitleSize(0.15);
  hist->GetXaxis()->SetTickLength(0.06);

  hist->GetYaxis()->SetLabelSize(0.15);
  hist->GetYaxis()->SetTitleOffset(0.7);
  hist->GetYaxis()->SetTitleSize(0.16);
  hist->GetYaxis()->SetTitle("#splitline{Normalised}{Residuals}");
  hist->GetYaxis()->SetNdivisions(204);

  hist->SetMaximum(+3.99);
  hist->SetMinimum(-3.99);

  for( unsigned int i=0; i<normres_sigma.size(); i++ )
    {
      normres_sigma.at(i).SetX1(hist->GetXaxis()->GetXmin());
      normres_sigma.at(i).SetX2(hist->GetXaxis()->GetXmax());
    }

  normres_sigma.at(0).SetY1(0);
  normres_sigma.at(0).SetY2(0);
  normres_sigma.at(1).SetY1(1);
  normres_sigma.at(1).SetY2(1);
  normres_sigma.at(2).SetY1(-1);
  normres_sigma.at(2).SetY2(-1);
  normres_sigma.at(3).SetY1(2);
  normres_sigma.at(3).SetY2(2);
  normres_sigma.at(4).SetY1(-2);
  normres_sigma.at(4).SetY2(-2);

  normres_sigma.at(1).SetLineStyle(kDotted);
  normres_sigma.at(1).SetLineColor(kRed);
  normres_sigma.at(2).SetLineStyle(kDotted);
  normres_sigma.at(2).SetLineColor(kRed);
  normres_sigma.at(3).SetLineStyle(kDashed);
  normres_sigma.at(3).SetLineColor(kBlue);
  normres_sigma.at(4).SetLineStyle(kDashed);
  normres_sigma.at(4).SetLineColor(kBlue);
}

// Calculate normalised residuals
void norm_res( const TH1* h_data, const TH1* h_pdf, TH1* h_norm_res )
{
  if( h_data->GetXaxis()->GetNbins() != h_pdf->GetXaxis()->GetNbins() )
    {
      std::cerr << "ERROR: Number of coefficients incompatible with order. Check code." << std::endl;
      exit(1);
    }

  for( int i=0; i<h_data->GetXaxis()->GetNbins(); i++ )
    {
      double norm_res;
      if( h_data->GetBinError(i+1) != 0.0 )
	norm_res =
	  (h_data->GetBinContent(i+1) - h_pdf->GetBinContent(i+1))/
	  h_data->GetBinError(i+1);
      else
	norm_res =
	  (h_data->GetBinContent(i+1) - h_pdf->GetBinContent(i+1))/
	  1.0;
      h_norm_res->SetBinContent(i+1, norm_res);
      h_norm_res->SetBinError(i+1, 1.0);
    }
}


#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
