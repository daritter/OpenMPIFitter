//
//  Ks_eff.cc v1.0
//  E. White (Univ. of Cincinnati)
//  Determine overall efficiency correction factor and uncertainty for Ks
//

#include <string>
#include <iostream>
#include <cmath>

#include "Ks_eff.h"

#include "belle.h"

#include "belleCLHEP/Vector/ThreeVector.h"
#include "belleCLHEP/Matrix/Vector.h"
#include "helix/Helix.h"
#include <stdio.h>


#include "panther/panther.h"
#include MDST_H
#include "belleutil/debugout.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif
  using std::cout;
  using std::endl;

// Initialize static member for counting bins and total num
//int ks_eff::bin_count[][4] = { {0, 0, 0},
                            //{0, 0, 0},
                            //{0, 0, 0} };
//int ks_eff::global_count = 0;

// from http://belle.kek.jp/secured/wiki/doku.php?id=software:systematics:ks_systematic
const double Ks_eff::table_eff[][4] = { {1.0000, 1.0000, 1.0000, 1.0000},
                                  {0.9808, 0.9716, 0.9727, 0.9863},
                                  {0.9423, 0.9657, 0.9984, 0.9938} };

const double Ks_eff::table_err[][4] = { {.0375, .0375, .0375, .0375},
                            {.0252, .0087, .0075, .0104},
                            {.0384, .0230, .0118, .0122} };

Ks_eff::Ks_eff():bin_count{{0,0,0,0},{0,0,0,0},{0,0,0,0}},global_count(0){}

void Ks_eff::add(double plab, double costheta)
{

  // Momentum and cos(theta) indices
  int p_index;
  int cos_index;

  // check for momentum index
  if ( plab <= 0.5 ) p_index =  0;
  if ( (0.5 < plab) && (plab <= 1.5) ) p_index =  1;
  if ( 1.5 < plab ) p_index = 2;

  // check for cos(theta) index
  if ( costheta <= -0.5 ) cos_index = 0;
  if ( (-0.5 < costheta) &&  (costheta <= 0.) ) cos_index = 1;
  if ( (0. < costheta) &&  (costheta <= 0.5) ) cos_index = 2;
  if ( 0.5 < costheta ) cos_index = 3;

  _p_i = p_index;
  _c_i = cos_index;

  // Increment global counter and bin counter
  ++global_count;
  ++bin_count[p_index][cos_index];
}

// Altrnate Constructor
//Ks_eff::Ks_eff(HepLorentzVector p_ks, HepDouble cos_th)
//{

  //// Momentum and cos(theta) indices
  //int p_index;
  //int cos_index;

  //// check for momentum index
  //if ( p_ks.rho() <= 0.5 ) p_index =  0;
  //if ( (0.5 < p_ks.rho()) && (p_ks.rho() <= 1.5) ) p_index =  1;
  //if ( 1.5 < p_ks.rho() ) p_index = 2;

  //// check for cos(theta) index
  //if ( cos_th <= -0.5 ) cos_index = 0;
  //if ( (-0.5 < cos_th) &&  (cos_th <= 0.) ) cos_index = 1;
  //if ( (0. < cos_th) &&  (cos_th <= 0.5) ) cos_index = 2;
  //if ( 0.5 < cos_th ) cos_index = 3;


  //_p_i = p_index;
  //_c_i = cos_index;

  //// Increment global counter and bin counter
  //++global_count;
  //++bin_count[p_index][cos_index];
//}

double Ks_eff::get_thisEff() {
  // return efficiency factor for just a single event
  return table_eff[_p_i][_c_i];

}

double Ks_eff::get_thisErr() {
  // return error for just a single event
  return table_err[_p_i][_c_i];

}

int Ks_eff::get_totalNum() {
  //returns the total number of events considered up to this point
  return global_count;
}

// call this at the very end
double Ks_eff::total_effFactor() {
  //Weighs each factor by the statistics recorded in bin_count[][]

  // Initialize
  double p_err[3] = {0,0,0};
  double p_wgt[3] = {0,0,0};
  // For efficiency factors
  double p_eff[3] = {0,0,0};

  //count how many in each row (sum cosTheta bins)
  int row_total[3] = {0,0,0};
  for (int p=0; p<3; ++p) {
    for (int theta=0; theta<4; ++theta) { row_total[p] +=  bin_count[p][theta]; }
  }

  // loop over all indices, add weighted factor N_i / (sigma_i)^2
  for (int p=0; p<3; ++p) {
  // for each row
  double cos_w = 0.; //sum of frac/sig^2 for each row
  double c_wgt = 0; // sum of frac^2 for each row
  double eff_w = 0.; //sum of eff*frac/sig^2 for each row
  double eff_wgt = 0.; //normalized sum of 1/sig^2
    for (int theta=0; theta<4; ++theta) {
      double frac_cos = 0.;
      // get number in this bin, divide by total to get frac
      if (row_total[p] > 0) frac_cos = double(bin_count[p][theta])/double(row_total[p]);
      // (if == 0 then frac_cos stays zero
      double err_sq = table_err[p][theta]*table_err[p][theta];
      // sum 1/sig^2 with weight (frac)
      cos_w += frac_cos*frac_cos/err_sq;
      eff_w += table_eff[p][theta]*frac_cos/err_sq;
      eff_wgt += frac_cos/err_sq;
      c_wgt += frac_cos*frac_cos;
    }
    // now get p weighted err
    if (cos_w > 0) p_err[p] = sqrt(c_wgt)/sqrt(cos_w); // total err for this p
    if (eff_w > 0) p_eff[p] = eff_w/eff_wgt;
    p_wgt[p] = double(row_total[p])/double(global_count); // frac in this p region
  }

  // fraction of events in each momentum region
  double fr1 = p_wgt[0];
  double fr2 = p_wgt[1];
  double fr3 = p_wgt[2];
  // safety net
  double term1 = 0; if (p_err[0] > 0) term1 = fr1/(p_err[0]*p_err[0]);
  double eff_term1 = 0; if (p_err[0] > 0) eff_term1 = p_eff[0]*fr1/(p_err[0]*p_err[0]);
  double term2 = 0; if (p_err[1] > 0) term2 = fr2/(p_err[1]*p_err[1]);
  double eff_term2 = 0; if (p_err[1] > 0) eff_term2 = p_eff[1]*fr2/(p_err[1]*p_err[1]);
  double term3 = 0; if (p_err[2] > 0) term3 = fr3/(p_err[2]*p_err[2]);
  double eff_term3 = 0; if (p_err[2] > 0) eff_term3 = p_eff[2]*fr3/(p_err[2]*p_err[2]);


  // weigh errors frac/sig^2
  double total_w = term1 + term2 + term3;

  double total_eff = 0.;
  if (total_w > 0) total_eff = (eff_term1+eff_term2+eff_term3)/total_w;

  return total_eff;

}


double Ks_eff::total_errFactor() {

  //Weighs each error by corresponding statistics in bin_count[][]

  // Initialize
  double p_err[3] = {0,0,0};
  double p_wgt[3] = {0,0,0};

  //count how many in each row (sum cosTheta bins)
  int row_total[3] = {0,0,0};
  for (int p=0; p<3; ++p) {
    for (int theta=0; theta<4; ++theta) { row_total[p] +=  bin_count[p][theta]; }
  }

  // loop over all indices, add weighted factor N_i / (sigma_i)^2
  for (int p=0; p<3; ++p) {
  // for each row
  double cos_w = 0.; //sum of frac/sig^2 for each row
  double c_wgt = 0; // sum of frac^2 for each row
    for (int theta=0; theta<4; ++theta) {
      double frac_cos = 0.;
      // get number in this bin, divide by total to get frac
      if (row_total[p] > 0) frac_cos = double(bin_count[p][theta])/double(row_total[p]);
      // (if == 0 then frac_cos stays zero
      double err_sq = table_err[p][theta]*table_err[p][theta];
      // sum 1/sig^2 with weight (frac)
      cos_w += frac_cos*frac_cos/err_sq;
      c_wgt += frac_cos*frac_cos;
    }
    // now get p weighted err
    if (cos_w > 0) p_err[p] = (c_wgt)/sqrt(cos_w); // total err for this p
    p_wgt[p] = double(row_total[p])/double(global_count); // frac in this p region
  }

  // safety net
  double term1 = 0; if (p_err[0] > 0) term1 = p_wgt[0]*p_wgt[0]/(p_err[0]*p_err[0]);
  double term2 = 0; if (p_err[1] > 0) term2 = p_wgt[1]*p_wgt[1]/(p_err[1]*p_err[1]);
  double term3 = 0; if (p_err[2] > 0) term3 = p_wgt[2]*p_wgt[2]/(p_err[2]*p_err[2]);

  // weigh errors frac/sig^2
  double total_w = term1 + term2 + term3;
  double total_err = 0.;
  if (total_w > 0) total_err = (p_wgt[0]*p_wgt[0] + p_wgt[1]*p_wgt[1] + p_wgt[2]*p_wgt[2])/sqrt(total_w);

  return total_err;

}

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

