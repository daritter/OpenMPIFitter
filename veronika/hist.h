#ifndef Hist_H_
#define Hist_H_

#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TLine.h"
#include "TCanvas.h"
#include "Data.h"

using namespace ROOT::Minuit2;

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

// Canvas style
void set_style();

// Histogram style
void hist_style( TH1* hist, const std::string& title_x,
		 const std::string& units="" );

// 1D Histogram PDF
double hist1d( const double& x,
	       const TH1D* h_ref );
double norm_hist1d( const double& x_ll, const double& x_ul,
		    const TH1D* h_ref );
void load_hist1d( std::ifstream& datafile, TH1D* h_ );
void symm_hist1d( TH1D* h_ );
double get_chisq( const TH1* h_data, const TH1* h_pdf,
		  const unsigned int& free_par );

// 2D Histogram PDF
double hist2d( const double& x, const double& y,
	       const TH2D* h_ref );
double norm_hist2d( const double& x_ll, const double& x_ul,
		    const double& y_ll, const double& y_ul,
		    const TH2D* h_ref );

// 3D Histogram PDF
double hist3d( const double& x, const double& y, const double& z,
	       const TH3D* h_ref );
double norm_hist3d( const double& x_ll, const double& x_ul,
		    const double& y_ll, const double& y_ul,
		    const double& z_ll, const double& z_ul,
		    const TH3D* h_ref );

// Load histograms and normalise
void load_dehist( std::vector<Data>& data, TH1D* h_de_svd1, TH1D* h_de_svd2 );
void load_ma1hist( std::vector<Data>& data,
		   TH1D* h_ma1_svd1, TH1D* h_ma1_svd2 );
void load_ha1hist( std::vector<Data>& data,
		   TH1D* h_ha1_svd1, TH1D* h_ha1_svd2 );
void load_dema1hist( std::vector<Data>& data,
		     TH2D* h_dema1_svd1, TH2D* h_dema1_svd2 );
void load_rbinhist( std::vector<Data>& data,
		    TH1D* h_rbin_svd1, TH1D* h_rbin_svd2 );
void normres_style( TCanvas* canvas );

void normres_style( TH1* hist, std::vector<TLine>& normres_sigma,
		    const std::string& title_x, const std::string& units="" );
void norm_res( const TH1* h_data, const TH1* h_pdf, TH1* h_norm_res );
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

#endif //Hist_H_
