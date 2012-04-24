#ifndef Wtag_H_
#define Wtag_H_

#include "belle.h"
#include "tatami/tatami.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TRandom3.h"

#include "Minuit2/FunctionMinimum.h"

#include "constant.h"
#include "Data.h"

using namespace ROOT::Minuit2;

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

//Wrong tag fraction
const unsigned int bins_rbin = 7;
const double rbin_bound[] = {0.0, 0.1, 0.25, 0.5, 0.625, 0.75, 0.875, 1.0};

//Load and normalise r-bin dependent histograms
void load_fdhist_rbin( std::vector<Data>& data,
		       std::vector<TH1D*>& h_fd_svd1,
		       std::vector<TH1D*>& h_fd_svd2,
		       const unsigned int& id );

const int set_rbin(const double r);

double set_wtag( const unsigned int& exp_no, const unsigned int& mc_type,
		 const int& rbin );
double set_dwtag( const unsigned int& exp_no, const unsigned int& mc_type,
		  const int& rbin );

const double set_wtag_uerr(const int expno, const int rbin);
const double set_wtag_lerr(const int expno, const int rbin);
const double set_delta_wtag_uerr(const int expno, const int rbin);
const double set_delta_wtag_lerr(const int expno, const int rbin);

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

#endif //Wtag_H_
