#ifndef Hist_H_
#define Hist_H_

#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"

// Canvas style
void set_style( TStyle* style );

// Histogram style
void hist_style( TH1* hist, const std::string& title_x );
void hist_style( TH2* hist, const std::string& title_x,
		 const std::string& title_y );

// Calculate chi-square
double get_chisq( const TH1* h_data, const TH1* h_pdf,
		  const unsigned int& free_par );
double get_chisq2D( const TH2* h_data, const TH2* h_pdf,
		  const unsigned int& free_par );
#endif //Hist_H_
