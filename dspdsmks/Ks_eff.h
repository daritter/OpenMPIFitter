#ifndef KS_EFF_H
#define KS_EFF_H

#include <cstring>
#include "belle.h"
#include <stdio.h>
#include "particle/utility.h"
#include "particle/combination.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

//class Particle;

class Ks_eff
{

public:
  Ks_eff();
  // Recommended cnstructor
  //Ks_eff(Particle& ks);
  //Ks_eff(HepLorentzVector p_ks, HepDouble cos_th); //alternate
  //virtual ~Ks_eff() {};
  void add(double plab, double costheta);
  double total_effFactor(); //total eff correction after all events
  double total_errFactor(); //total uncertainty after all events
  double get_thisEff(); //returns efficiency factor for this single event
  double get_thisErr(); //returns error for just this event
  int get_totalNum(); // return the total number of events up to this point
                       // (total number being used in weighing the eff. and errors)

private:
  int bin_count[3][4];// counts events in momentum/cosTheta bins
  int global_count; // keeps track of total num of events (sum of all bin_count[][4])
  // index values
  int _p_i; //
  int _c_i;
  //void _calc_errs();
  //  static int bin_count[3];
  static const double table_err[][4];
  static const double table_eff[][4];

};

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

#endif
