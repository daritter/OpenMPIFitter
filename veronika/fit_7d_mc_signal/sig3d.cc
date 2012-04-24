#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnPrint.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH3.h"
#include "TLatex.h"
#include "belle.h"
#include "tatami/tatami.h"
#include "TLegend.h"

#include "../wtag.h"
#include "../hist.h"
#include "../func.h"
#include "../cheb.h"
#include "../tools.h"
#include "Sig3DFcn.h"

#include "MPIFitter.h"

using namespace Belle;
using namespace Belle::ROOT::Minuit2;

void plot_histogram( Sig3DFcn& fcn, const std::vector<double>& par,
             TCanvas *canvas1, TCanvas *canvas2, const int& rbin, std::vector<Data> data, const std::string& mode );

struct FitRoutine {
    template<class FCN> int operator()(FCN& fcn) {
        //Fitting here

  int n_par_wks = 106;
  int n_par_wkp = 106;
  // Load parameters
  MnUserParameters mn_param;

  std::ifstream parfile;
  parfile.open("veronika/fit_7d_mc_signal/initial_parameters.dat");
  load_par( parfile, mn_param );

     mn_param.Fix(10);
     mn_param.Fix(10+n_par_wks);
     mn_param.Fix(10+2*n_par_wks);
     mn_param.Fix(10+2*n_par_wks+n_par_wkp);
     mn_param.Fix(11);
     mn_param.Fix(11+n_par_wks);
     mn_param.Fix(11+2*n_par_wks);
     mn_param.Fix(11+2*n_par_wks+n_par_wkp);

  /* for(unsigned int i = 0; i <10; i++){
     mn_param.Fix(i+2*n_par_wks);
     mn_param.Fix(i+2*n_par_wks+n_par_wkp);
   }    */
    for(unsigned int i = 0; i <106; i++){
     mn_param.Fix(i);
     mn_param.Fix(i+n_par_wks);
     mn_param.Fix(i+2*n_par_wks);
     mn_param.Fix(i+2*n_par_wks+n_par_wkp);
   }
// Delta E
  for(unsigned int i = 0; i <10; i++){
     mn_param.Release(i);
     mn_param.Release(i+n_par_wks);
   }
   for(unsigned int i = 0; i <2; i++){
     mn_param.Release(i+2*n_par_wks);
     mn_param.Release(i+2*n_par_wks+n_par_wkp);
   }

// M omega
   for(unsigned int i = 68; i <78; i++){
    mn_param.Release(i);
    mn_param.Release(i+n_par_wks);
   }
   for(unsigned int i = 76; i <78; i++){
    mn_param.Release(i+2*n_par_wks);
    mn_param.Release(i+2*n_par_wks+n_par_wkp);
   }

// H omega
    mn_param.Release(79);
    mn_param.Release(79+n_par_wks);
    mn_param.Release(81);
    mn_param.Release(81+n_par_wks);

// Mbc
   for(unsigned int i = 82; i <90; i++){
    mn_param.Release(i);
    mn_param.Release(i+n_par_wks);
   }
   for(unsigned int i = 82; i <84; i++){
    mn_param.Release(i+2*n_par_wks);
    mn_param.Release(i+2*n_par_wks+n_par_wkp);
   }

// Fisher

   for(unsigned int i = 12; i <68; i++){
     mn_param.Release(i);
     mn_param.Release(i+n_par_wks);
   }
   for(unsigned int i = 12; i <14; i++){
     mn_param.Release(i+2*n_par_wks);
     mn_param.Release(i+2*n_par_wks+n_par_wkp);
   }
   for(unsigned int i = 20; i <22; i++){
     mn_param.Release(i+2*n_par_wks);
     mn_param.Release(i+2*n_par_wks+n_par_wkp);
   }
  for(unsigned int i = 28; i <30; i++){
     mn_param.Release(i+2*n_par_wks);
     mn_param.Release(i+2*n_par_wks+n_par_wkp);
   }
   for(unsigned int i = 36; i <38; i++){
     mn_param.Release(i+2*n_par_wks);
     mn_param.Release(i+2*n_par_wks+n_par_wkp);
   }
   for(unsigned int i = 44; i <46; i++){
     mn_param.Release(i+2*n_par_wks);
     mn_param.Release(i+2*n_par_wks+n_par_wkp);
   }
  for(unsigned int i = 52; i <54; i++){
     mn_param.Release(i+2*n_par_wks);
     mn_param.Release(i+2*n_par_wks+n_par_wkp);
   }
    for(unsigned int i = 60; i <62; i++){
     mn_param.Release(i+2*n_par_wks);
     mn_param.Release(i+2*n_par_wks+n_par_wkp);
   }
/**/
   for(unsigned int i = 93; i <99; i++){
     mn_param.Release(i);
     mn_param.SetLimits(i, 0., 1.);
     mn_param.Release(i+n_par_wks);
     mn_param.Release(i+2*n_par_wks);
     mn_param.SetLimits(i+2*n_par_wks, 0., 1.);
     mn_param.Release(i+2*n_par_wks+n_par_wkp);
     mn_param.SetLimits(i+2*n_par_wks+n_par_wkp, 0., 1.);
   }
     mn_param.Release(105);
     mn_param.Release(105+n_par_wks);
     mn_param.Release(105+2*n_par_wks);
     mn_param.Release(105+2*n_par_wks+n_par_wkp);


  for( unsigned int i=0; i<bins_rbin; i++ )
    {
      mn_param.SetLimits(12+8*i+0, h_fd_ll, h_fd_ul);
      mn_param.SetLimits(12+n_par_wks+8*i+0,        h_fd_ll, h_fd_ul);
      mn_param.SetLimits(12+2*n_par_wks+8*i+0,        h_fd_ll, h_fd_ul);
      mn_param.SetLimits(12+2*n_par_wks+n_par_wkp+8*i+0,        h_fd_ll, h_fd_ul);
      mn_param.SetLimits(12+8*i+1, 1.0, 10.);
      mn_param.SetLimits(12+2*n_par_wks+8*i+1, 1.0, 10.);
      mn_param.SetLimits(12+2*n_par_wks+n_par_wkp+8*i+1, 1.0, 10.);
      //mn_param.Fix(12+8*i+2);
      //mn_param.Fix(12+8*i+3);
      mn_param.SetLimits(12+8*i+2, h_fd_ll, h_fd_ul);
      mn_param.SetLimits(12+8*i+3, 1.0, 10.);
      mn_param.SetLimits(12+8*i+4, h_fd_ll, h_fd_ul);
      mn_param.SetLimits(12+n_par_wks+8*i+4,        h_fd_ll, h_fd_ul);
      mn_param.SetLimits(12+8*i+5, 1.0, 10.);
      mn_param.SetLimits(12+8*i+6, 0.0, 1.);
      mn_param.SetLimits(12+n_par_wks+8*i+6, 0.0, 1.);
     // mn_param.Fix(12+8*i+7);
      mn_param.SetLimits(12+8*i+7, 0.0, 1.);
      mn_param.SetLimits(12+n_par_wks+8*i+7, 0.0, 1.);

    }
     mn_param.SetLimits("1sigwksfd_s16", 0., 3.);
     mn_param.SetLimits("1sigwkpfd_s16", 0., 3.);
     mn_param.SetLimits("1sigwksfd_s26", 0., 3.);
     mn_param.SetLimits("1sigwksfd_s36", 0., 3.);
     mn_param.SetLimits("2sigwkpfd_s16", 0., 3.);
     mn_param.SetLimits("1sigwkswm_s3", 0., 4.);
     mn_param.SetLimits("2sigwkswm_m3", -0.5, 0.5);
     mn_param.SetLimits("2sigwkswm_m2", -0.5, 0.5);
     mn_param.SetLimits("2sigwkswm_s1", 0., 4.);
     mn_param.SetLimits("2sigwkswm_s2", 0., 4.);
     mn_param.SetLimits("1sigwkswm_corr_s", 0., 1.);
     mn_param.SetLimits("2sigwkswm_corr_s", 0., 1.);
     mn_param.SetLimits("1sigwkpwm_corr_s", 0., 1.);
     mn_param.SetLimits("2sigwkpwm_corr_s", 0., 1.);

//   //const double h_mbc_ll = 5.259, h_mbc_ul = 5.291;

   //mn_param.Release(73);
  std::cout << "Initial parameter state:" << mn_param << std::endl;


  // Minimisation
  MnMigrad migrad( fcn, mn_param, 2 ); //Strategy 2
  FunctionMinimum min = migrad();
//   min = migrad();
//   min = migrad();
//   min = migrad();
//   min = migrad();

  // Print the result
  std::cout << min << std::endl;
  std::cout << "Function Minimum: " << std::setprecision(10)
	    << min.Fval() << std::endl;

  //if( min.IsValid() == false )
  //exit(1);params
  std::vector<double> par = min.UserParameters().Params();

  //std::vector<double> par = mn_param.Params();

  // Save file
  std::ofstream savefile;
  savefile.open("fit_result_svd.par", std::ios::out);
  for( unsigned int i=0; i<mn_param.Params().size(); i++ )
    savefile << setiosflags(std::ios::left) << std::setprecision(10)
             << std::setw(20) << mn_param.Name(i)
             << std::setw(20) << min.UserParameters().Params().at(i)
             << std::setw(20) << min.UserParameters().Errors().at(i)
             << std::endl;
  savefile.close();

  return 0;

    }
};

int main( int argc, char* argv[] )
{
  // Load data
    std::cout << "Aye, loading thy data" << std::endl;
  std::ifstream datafile11, datafile12, datafile21, datafile22;
  std::string path = "/remote/pcbelle07/veronika/belle/fit_combined";
  datafile11.open((path+"/svd_mc_signal_wks.txt").c_str());
  datafile12.open((path+"/svd_mc_signal_wks_thrust.txt").c_str());
  datafile21.open((path+"/svd_mc_signal_wkp.txt").c_str());
  datafile22.open((path+"/svd_mc_signal_wkp_thrust.txt").c_str());

  std::vector<Data> data10;
  FillPlot(datafile11, datafile12, data10, 1, 0, 1);

  std::vector<Data> data11;
  FillPlot(datafile21, datafile22, data11, 1, 0, 1);

 std::cout << "Data fully loaded capt'n" << std::endl;

  // Initialize your Fcn
  Sig3DFcn fcn( data10, data11 );
  // Get a fitting routine
  FitRoutine fitter;
  // Call the MPI Fitting core and return the result
  MPIFitter core;
  return core.run(fitter, fcn);

}

