#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>

#include "hist.h"

// Canvas style
void set_style( TStyle* style )
{
  gStyle->SetPadLeftMargin(0.22);
  gStyle->SetPadBottomMargin(0.20);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetHistLineWidth(2);
}

// Histogram style
void hist_style( TH1* hist, const std::string& title_x )
{
  hist->GetXaxis()->SetLabelSize(0.07);
  hist->GetXaxis()->SetLabelOffset(0.03);
  hist->GetXaxis()->SetTitleOffset(1.1);
  hist->GetXaxis()->SetTitleSize(0.08);
  hist->GetXaxis()->SetNdivisions(50205, kFALSE);
  hist->GetYaxis()->SetLabelSize(0.07);
  hist->GetYaxis()->SetTitleOffset(1.4);
  hist->GetYaxis()->SetTitleSize(0.08);

  hist->GetXaxis()->SetTitle(title_x.c_str());

  const double bin_width =
    (hist->GetXaxis()->GetXmax() - hist->GetXaxis()->GetXmin())/
    static_cast<double>(hist->GetXaxis()->GetNbins());
  std::stringstream ss_ytitle;
  ss_ytitle << bin_width;
  std::string s_ytitle = "Events / (" + ss_ytitle.str() + ")";
  hist->GetYaxis()->SetTitle(s_ytitle.c_str());
}

void hist_style( TH2* hist, const std::string& title_x,
		 const std::string& title_y )
{
  hist->GetXaxis()->SetLabelSize(0.07);
  hist->GetXaxis()->SetLabelOffset(0.03);
  hist->GetXaxis()->SetTitleOffset(1.1);
  hist->GetXaxis()->SetTitleSize(0.08);
  hist->GetXaxis()->SetNdivisions(50205, kFALSE);
  hist->GetYaxis()->SetLabelSize(0.07);
  hist->GetYaxis()->SetTitleOffset(1.4);
  hist->GetYaxis()->SetTitleSize(0.08);
  hist->GetYaxis()->SetNdivisions(50205, kFALSE);

  hist->GetXaxis()->SetTitle(title_x.c_str());
  hist->GetYaxis()->SetTitle(title_y.c_str());
}

// Calculate chi-square
double get_chisq( const TH1* h_data, const TH1* h_pdf,
		  const unsigned int& free_par )
{
  if( h_data->GetXaxis()->GetNbins() != h_pdf->GetXaxis()->GetNbins() )
    {
      std::cout << "Error: Data and PDF histograms must have equal bins"
		<< std::endl;
      std::exit(1);
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

double get_chisq2D( const TH2* h_data, const TH2* h_pdf,
		  const unsigned int& free_par )
{
  if( h_data->GetXaxis()->GetNbins() != h_pdf->GetXaxis()->GetNbins() )
    {
      std::cout << "Error: Data and PDF histograms must have equal bins"
		<< std::endl;
      std::exit(1);
    }

  double chisq = 0.0;
  unsigned int bins = 0;
  for( int i=0; i<h_data->GetXaxis()->GetNbins(); i++ )
    {
	  for( int j=0; j<h_data->GetYaxis()->GetNbins(); j++ )
   	{
      const double data = h_data->GetBinContent(i+1, j+1);
      const double pdf  = h_pdf->GetBinContent(i+1, j+1);
      const double err  = sqrt(data);

      	if( data != 0.0 )
		{
	  chisq += (data-pdf)*(data-pdf)/err/err;
	  bins++;
		}
  	  }
	}

  return chisq/static_cast<double>(bins-free_par);
}
