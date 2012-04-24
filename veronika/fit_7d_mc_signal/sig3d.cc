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

using namespace Belle;
using namespace Belle::ROOT::Minuit2;

void plot_histogram( Sig3DFcn& fcn, const std::vector<double>& par,
             TCanvas *canvas1, TCanvas *canvas2, const int& rbin, std::vector<Data> data, const std::string& mode );

int main( int argc, char* argv[] )
{
  std::vector<TLine> normres_sigma(5);
  set_style();

  int n_par_wks = 106;
  int n_par_wkp = 106;
  // Load data
  std::ifstream datafile11, datafile12, datafile21, datafile22;
  datafile11.open("../svd_mc_signal_wks.txt");
  datafile12.open("../svd_mc_signal_wks_thrust.txt");
  datafile21.open("../svd_mc_signal_wkp.txt");
  datafile22.open("../svd_mc_signal_wkp_thrust.txt");
  
  std::vector<Data> data10;
  FillPlot(datafile11, datafile12, data10, 1, 0, 1);

  std::vector<Data> data11;
  FillPlot(datafile21, datafile22, data11, 1, 0, 1);  
  // Load parameters
  MnUserParameters mn_param;

  std::ifstream parfile;
  parfile.open("initial_parameters.dat");
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

  // Initialize your Fcn
  Sig3DFcn fcn( data10, data11 );

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

  const unsigned int bins_oh = 50, bins_dt = 15;
 // const double h_mbc_ll = 5.265, h_mbc_ul = 5.29;
 /* 
  //Plot de  
  TH1D *h_de_svd1 = new TH1D( "h_de_svd1", "", bins_de, h_de_ll, h_de_ul);
  TH1D *h_de_svd2 = new TH1D( "h_de_svd2", "", bins_de, h_de_ll, h_de_ul);
  for( unsigned int i=0; i<data10.size(); i++ ){
    if( data10[i].GetExpNo() < 29 )
      h_de_svd1->Fill( data10[i].GetDE());
    else
      h_de_svd2->Fill( data10[i].GetDE());
  }
  
  TH1D *h_de_svd1_pdf = new TH1D( "h_de_svd1_pdf", "", bins_de, h_de_ll, h_de_ul );    
  TH1D *h_de_svd2_pdf = new TH1D( "h_de_svd2_pdf", "", bins_de, h_de_ll, h_de_ul );
   
   const double de_width =
    ( h_de_ul - h_de_ll ) / static_cast<double>(bins_de);
    
   for( unsigned int k=0; k < bins_de; k++ )
      {
     	 const double h_de_svd1 =
		h_de_ll + ((0.5+static_cast<double>(k))*de_width);     	 
	 const double h_de_svd2 =
		h_de_ll + ((0.5+static_cast<double>(k))*de_width);
     	 double h_pdf_de_svd1 = 0.0;
     	 double h_pdf_de_svd2 = 0.0;
     
	for( unsigned int l=0; l<data10.size(); l++ ){	
	      if( data10[l].GetExpNo() < 29 )
		h_pdf_de_svd1 += fcn.wksDE(data10[l], h_de_svd1, par)/fcn.NormwksDE("svd1", par);
	      else
		h_pdf_de_svd2 += fcn.wksDE(data10[l], h_de_svd2, par)/fcn.NormwksDE("svd2", par);      
	 }
      h_de_svd1_pdf->Fill(h_de_svd1, h_pdf_de_svd1);
      h_de_svd2_pdf->Fill(h_de_svd2, h_pdf_de_svd2);
      }
 
  h_de_svd1_pdf->Scale(de_width);
  h_de_svd2_pdf->Scale(de_width);
 
  
    //Plot de  
 TH1D *h_de_svd1_wkp = new TH1D( "h_de_svd1_wkp", "", bins_de, h_de_ll, h_de_ul);
  TH1D *h_de_svd2_wkp = new TH1D( "h_de_svd2_wkp", "", bins_de, h_de_ll, h_de_ul);
  for( unsigned int i=0; i<data11.size(); i++ ){
    if( data11[i].GetExpNo() < 29 )
      h_de_svd1_wkp->Fill( data11[i].GetDE());
    else
      h_de_svd2_wkp->Fill( data11[i].GetDE());
  }
  
  TH1D *h_de_svd1_pdf_wkp = new TH1D( "h_de_svd1_pdf_wkp", "", bins_de, h_de_ll, h_de_ul );    
  TH1D *h_de_svd2_pdf_wkp = new TH1D( "h_de_svd2_pdf_wkp", "", bins_de, h_de_ll, h_de_ul );

    
   for( unsigned int k=0; k < bins_de; k++ )
      {
     	 const double h_de_svd1_wkp =
		h_de_ll + ((0.5+static_cast<double>(k))*de_width);     	 
	 const double h_de_svd2_wkp =
		h_de_ll + ((0.5+static_cast<double>(k))*de_width);
     	 double h_pdf_de_svd1_wkp = 0.0;
     	 double h_pdf_de_svd2_wkp = 0.0;
     
	for( unsigned int l=0; l<data11.size(); l++ ){	
	      if( data11[l].GetExpNo() < 29 )
		h_pdf_de_svd1_wkp += fcn.wkpDE(data11[l], h_de_svd1_wkp, par)/fcn.NormwkpDE("svd1", par);
	      else
		h_pdf_de_svd2_wkp += fcn.wkpDE(data11[l], h_de_svd2_wkp, par)/fcn.NormwkpDE("svd2", par);      
	 }
      h_de_svd1_pdf_wkp->Fill(h_de_svd1_wkp, h_pdf_de_svd1_wkp);
      h_de_svd2_pdf_wkp->Fill(h_de_svd2_wkp, h_pdf_de_svd2_wkp);
      }
 
  h_de_svd1_pdf_wkp->Scale(de_width);
  h_de_svd2_pdf_wkp->Scale(de_width);
 
  
  //Plot om
  TH1D *h_om_svd1 = new TH1D( "h_om_svd1", "", bins_om, h_momega_ll, h_momega_ul);
  TH1D *h_om_svd2 = new TH1D( "h_om_svd2", "", bins_om, h_momega_ll, h_momega_ul);
  for( unsigned int i=0; i<data10.size(); i++ ){
    if( data10[i].GetExpNo() < 29 )
      h_om_svd1->Fill( data10[i].GetMomega());
    else
      h_om_svd2->Fill( data10[i].GetMomega());
  }
  
  TH1D *h_om_svd1_pdf = new TH1D( "h_om_svd1_pdf", "", bins_om, h_momega_ll, h_momega_ul );    
  TH1D *h_om_svd2_pdf = new TH1D( "h_om_svd2_pdf", "", bins_om, h_momega_ll, h_momega_ul );
   
   const double om_width =
    ( h_momega_ul - h_momega_ll ) / static_cast<double>(bins_om);
    
   for( unsigned int k=0; k < bins_om; k++ )
      {
     	 const double h_om_svd1 =
		h_momega_ll + ((0.5+static_cast<double>(k))*om_width);     	 
	 const double h_om_svd2 =
		h_momega_ll + ((0.5+static_cast<double>(k))*om_width);
     	 double h_pdf_om_svd1 = 0.0;
     	 double h_pdf_om_svd2 = 0.0;
     
	for( unsigned int l=0; l<data10.size(); l++ ){	
	      if( data10[l].GetExpNo() < 29 )
		h_pdf_om_svd1 += fcn.wksMw(data10[l], h_om_svd1, par)/fcn.NormwksMw("svd1", par);
	      else
		h_pdf_om_svd2 += fcn.wksMw(data10[l], h_om_svd2, par)/fcn.NormwksMw("svd2", par);
	 }
      h_om_svd1_pdf->Fill(h_om_svd1, h_pdf_om_svd1);
      h_om_svd2_pdf->Fill(h_om_svd2, h_pdf_om_svd2);
      }
 
  h_om_svd1_pdf->Scale(om_width);
  h_om_svd2_pdf->Scale(om_width);  
  
  TH1D *h_om_svd1_wkp = new TH1D( "h_om_svd1_wkp", "", bins_om, h_momega_ll, h_momega_ul);
  TH1D *h_om_svd2_wkp = new TH1D( "h_om_svd2_wkp", "", bins_om, h_momega_ll, h_momega_ul);
  for( unsigned int i=0; i<data11.size(); i++ ){
    if( data11[i].GetExpNo() < 29 )
      h_om_svd1_wkp->Fill( data11[i].GetMomega());
    else
      h_om_svd2_wkp->Fill( data11[i].GetMomega());
  }
  
  TH1D *h_om_svd1_pdf_wkp = new TH1D( "h_om_svd1_pdf_wkp", "", bins_om, h_momega_ll, h_momega_ul );    
  TH1D *h_om_svd2_pdf_wkp = new TH1D( "h_om_svd2_pdf_wkp", "", bins_om, h_momega_ll, h_momega_ul );
    
   for( unsigned int k=0; k < bins_om; k++ )
      {
     	 const double h_om_svd1_wkp =
		h_momega_ll + ((0.5+static_cast<double>(k))*om_width);     	 
	 const double h_om_svd2_wkp =
		h_momega_ll + ((0.5+static_cast<double>(k))*om_width);
     	 double h_pdf_om_svd1_wkp = 0.0;
     	 double h_pdf_om_svd2_wkp = 0.0;
     
	for( unsigned int l=0; l<data11.size(); l++ ){	
	      if( data11[l].GetExpNo() < 29 )
		h_pdf_om_svd1_wkp += fcn.wksMw(data11[l], h_om_svd1_wkp, par)/fcn.NormwksMw("svd1", par);
	      else
		h_pdf_om_svd2_wkp += fcn.wksMw(data11[l], h_om_svd2_wkp, par)/fcn.NormwksMw("svd2", par);
	 }
      h_om_svd1_pdf_wkp->Fill(h_om_svd1_wkp, h_pdf_om_svd1_wkp);
      h_om_svd2_pdf_wkp->Fill(h_om_svd2_wkp, h_pdf_om_svd2_wkp);
      }
 
  h_om_svd1_pdf_wkp->Scale(om_width);
  h_om_svd2_pdf_wkp->Scale(om_width);  
  */

  // Plot de mw
  TH2D *p_de_mw = new TH2D( "p_de_mw",  "", bins_de,  h_de_ll,  h_de_ul, bins_momega, h_momega_ll, h_momega_ul  );  
  TH2D *p_de_mw_sig = new TH2D( "p_de_mw_sig",  "", bins_de,  h_de_ll,  h_de_ul, bins_momega, h_momega_ll, h_momega_ul  );  

  TH1D *h_de_sig = new TH1D( "h_de_sig", "", bins_de, h_de_ll, h_de_ul);
  TH1D *h_om_sig = new TH1D( "h_om_sig", "", bins_momega, h_momega_ll, h_momega_ul);

  for( unsigned int i=0; i<data10.size(); i++ ){
      h_de_sig->Fill( data10[i].GetDE());
      h_om_sig->Fill( data10[i].GetMomega());
      p_de_mw->Fill(data10[i].GetDE(), data10[i].GetMomega());
   }
 
  TH1D *h_de_pdf = new TH1D( "h_de_pdf", "", bins_de, h_de_ll, h_de_ul );     
  TH1D *h_om_pdf = new TH1D( "h_om_pdf", "", bins_momega, h_momega_ll, h_momega_ul );     
 
  
	for( unsigned int k=0; k < bins_de; k++ ){
	  const double h_de = h_de_sig->GetXaxis()->GetBinCenter(k+1);
	  for( unsigned int m=0; m < bins_momega; m++ ){
	      const double h_om = h_om_sig->GetXaxis()->GetBinCenter(m+1);
	      
	      double h_sig_pdf = 0.0;
              fcn.PlotMwPDF( h_de, h_om, h_sig_pdf, par);
	      h_de_pdf->Fill(h_de, h_sig_pdf);
	      h_om_pdf->Fill(h_om, h_sig_pdf);
	      p_de_mw_sig->Fill(h_de, h_om, h_sig_pdf);
	  }
	}
	
      const double sf_mw = h_de_pdf->GetBinWidth(0)*h_om_pdf->GetBinWidth(0);  
      h_de_pdf->Scale(sf_mw);
      h_om_pdf->Scale(sf_mw);
      
    TH1D *h_de_mw_proj1 = p_de_mw->ProjectionY("de_mw_py1", 1., 10.);   
    TH1D *h_de_mw_proj2 = p_de_mw->ProjectionY("de_mw_py2", 11., 20.);   
    TH1D *h_de_mw_proj3 = p_de_mw->ProjectionY("de_mw_py3", 21., 30.);
    TH1D *h_de_mw_proj4 = p_de_mw->ProjectionY("de_mw_py4", 31., 40.);   
    TH1D *h_de_mw_proj5 = p_de_mw->ProjectionY("de_mw_py5", 41., 50.);
    TH1D *p_de_mw_proj1 = p_de_mw_sig->ProjectionY("de_mw_py1_sig", 1., 10.);   
    TH1D *p_de_mw_proj2 = p_de_mw_sig->ProjectionY("de_mw_py2_sig", 11., 20.);   
    TH1D *p_de_mw_proj3 = p_de_mw_sig->ProjectionY("de_mw_py3_sig", 21., 30.);
    TH1D *p_de_mw_proj4 = p_de_mw_sig->ProjectionY("de_mw_py4_sig", 31., 40.);   
    TH1D *p_de_mw_proj5 = p_de_mw_sig->ProjectionY("de_mw_py5_sig", 41., 50.);
    TH1D *p_de_mw_proj1x = p_de_mw->ProjectionX("de_mw_px1_sig", 1., 10.);   
    TH1D *p_de_mw_proj2x = p_de_mw->ProjectionX("de_mw_px2_sig", 11., 20.);   
    TH1D *p_de_mw_proj3x = p_de_mw->ProjectionX("de_mw_px3_sig", 21., 30.);
    TH1D *p_de_mw_proj4x = p_de_mw->ProjectionX("de_mw_px4_sig", 31., 40.);   
    TH1D *p_de_mw_proj5x = p_de_mw->ProjectionX("de_mw_px5_sig", 41., 50.);
    p_de_mw_proj1->Scale(p_de_mw_proj1->GetBinWidth(0)*p_de_mw_proj1x->GetBinWidth(0));
    p_de_mw_proj2->Scale(p_de_mw_proj2->GetBinWidth(0)*p_de_mw_proj2x->GetBinWidth(0));
    p_de_mw_proj3->Scale(p_de_mw_proj3->GetBinWidth(0)*p_de_mw_proj3x->GetBinWidth(0));
    p_de_mw_proj4->Scale(p_de_mw_proj4->GetBinWidth(0)*p_de_mw_proj4x->GetBinWidth(0));
    p_de_mw_proj5->Scale(p_de_mw_proj5->GetBinWidth(0)*p_de_mw_proj5x->GetBinWidth(0));    

    TH1D *r_de1 = new TH1D( "r_de1", "", bins_momega, h_momega_ll, h_momega_ul );
    norm_res( h_de_mw_proj1, p_de_mw_proj1, r_de1 );
    TH1D *r_de2 = new TH1D( "r_de2", "", bins_momega, h_momega_ll, h_momega_ul );
    norm_res( h_de_mw_proj2, p_de_mw_proj2, r_de2 );
    TH1D *r_de3 = new TH1D( "r_de3", "", bins_momega, h_momega_ll, h_momega_ul );
    norm_res( h_de_mw_proj3, p_de_mw_proj3, r_de3 );
    TH1D *r_de4 = new TH1D( "r_de4", "", bins_momega, h_momega_ll, h_momega_ul );
    norm_res( h_de_mw_proj4, p_de_mw_proj4, r_de4 );
    TH1D *r_de5 = new TH1D( "r_de5", "", bins_momega, h_momega_ll, h_momega_ul );
    norm_res( h_de_mw_proj5, p_de_mw_proj5, r_de5 );

    
  TH2D *p_de_mw_wkp = new TH2D( "p_de_mw_wkp",  "", bins_de,  h_de_ll,  h_de_ul, bins_momega, h_momega_ll, h_momega_ul  );  
  TH2D *p_de_mw_sig_wkp = new TH2D( "p_de_mw_sig_wkp",  "", bins_de,  h_de_ll,  h_de_ul, bins_momega, h_momega_ll, h_momega_ul  );  

  TH1D *h_de_sig_wkp = new TH1D( "h_de_sig_wkp", "", bins_de, h_de_ll, h_de_ul);
  TH1D *h_om_sig_wkp = new TH1D( "h_om_sig_wkp", "", bins_momega, h_momega_ll, h_momega_ul);

  for( unsigned int i=0; i<data11.size(); i++ ){
      h_de_sig_wkp->Fill( data11[i].GetDE());
      h_om_sig_wkp->Fill( data11[i].GetMomega());
      p_de_mw_wkp->Fill(data11[i].GetDE(), data11[i].GetMomega());
   }
 
  TH1D *h_de_pdf_wkp = new TH1D( "h_de_pdf_wkp", "", bins_de, h_de_ll, h_de_ul );     
  TH1D *h_om_pdf_wkp = new TH1D( "h_om_pdf_wkp", "", bins_momega, h_momega_ll, h_momega_ul );     
 
  
	for( unsigned int k=0; k < bins_de; k++ ){
	  const double h_de = h_de_sig_wkp->GetXaxis()->GetBinCenter(k+1);
	  for( unsigned int m=0; m < bins_momega; m++ ){
	      const double h_om = h_om_sig_wkp->GetXaxis()->GetBinCenter(m+1);
	      
	      double h_sig_pdf_wkp = 0.0;
              fcn.PlotMwPDF_wkp( h_de, h_om, h_sig_pdf_wkp, par);
	      h_de_pdf_wkp->Fill(h_de, h_sig_pdf_wkp);
	      h_om_pdf_wkp->Fill(h_om, h_sig_pdf_wkp);
	      p_de_mw_sig_wkp->Fill(h_de, h_om, h_sig_pdf_wkp);
	  }
	}
	
      h_de_pdf_wkp->Scale(sf_mw);
      h_om_pdf_wkp->Scale(sf_mw);
      
    TH1D *h_de_mw_proj1_wkp = p_de_mw_wkp->ProjectionY("de_mw_py1_wkp", 1., 10.);   
    TH1D *h_de_mw_proj2_wkp = p_de_mw_wkp->ProjectionY("de_mw_py2_wkp", 11., 20.);   
    TH1D *h_de_mw_proj3_wkp = p_de_mw_wkp->ProjectionY("de_mw_py3_wkp", 21., 30.);
    TH1D *h_de_mw_proj4_wkp = p_de_mw_wkp->ProjectionY("de_mw_py4_wkp", 31., 40.);   
    TH1D *h_de_mw_proj5_wkp = p_de_mw_wkp->ProjectionY("de_mw_py5_wkp", 41., 50.);
    TH1D *p_de_mw_proj1_wkp = p_de_mw_sig_wkp->ProjectionY("de_mw_py1_sig_wkp", 1., 10.);   
    TH1D *p_de_mw_proj2_wkp = p_de_mw_sig_wkp->ProjectionY("de_mw_py2_sig_wkp", 11., 20.);   
    TH1D *p_de_mw_proj3_wkp = p_de_mw_sig_wkp->ProjectionY("de_mw_py3_sig_wkp", 21., 30.);
    TH1D *p_de_mw_proj4_wkp = p_de_mw_sig_wkp->ProjectionY("de_mw_py4_sig_wkp", 31., 40.);   
    TH1D *p_de_mw_proj5_wkp = p_de_mw_sig_wkp->ProjectionY("de_mw_py5_sig_wkp", 41., 50.);
    TH1D *p_de_mw_proj1x_wkp = p_de_mw_wkp->ProjectionX("de_mw_px1_sig_wkp", 1., 10.);   
    TH1D *p_de_mw_proj2x_wkp = p_de_mw_wkp->ProjectionX("de_mw_px2_sig_wkp", 11., 20.);   
    TH1D *p_de_mw_proj3x_wkp = p_de_mw_wkp->ProjectionX("de_mw_px3_sig_wkp", 21., 30.);
    TH1D *p_de_mw_proj4x_wkp = p_de_mw_wkp->ProjectionX("de_mw_px4_sig_wkp", 31., 40.);   
    TH1D *p_de_mw_proj5x_wkp = p_de_mw_wkp->ProjectionX("de_mw_px5_sig_wkp", 41., 50.);
    p_de_mw_proj1_wkp->Scale(p_de_mw_proj1_wkp->GetBinWidth(0)*p_de_mw_proj1x_wkp->GetBinWidth(0));
    p_de_mw_proj2_wkp->Scale(p_de_mw_proj2_wkp->GetBinWidth(0)*p_de_mw_proj2x_wkp->GetBinWidth(0));
    p_de_mw_proj3_wkp->Scale(p_de_mw_proj3_wkp->GetBinWidth(0)*p_de_mw_proj3x_wkp->GetBinWidth(0));
    p_de_mw_proj4_wkp->Scale(p_de_mw_proj4_wkp->GetBinWidth(0)*p_de_mw_proj4x_wkp->GetBinWidth(0));
    p_de_mw_proj5_wkp->Scale(p_de_mw_proj5_wkp->GetBinWidth(0)*p_de_mw_proj5x_wkp->GetBinWidth(0));    

    TH1D *r_de1_wkp = new TH1D( "r_de1_wkp", "", bins_momega, h_momega_ll, h_momega_ul );
    norm_res( h_de_mw_proj1_wkp, p_de_mw_proj1_wkp, r_de1_wkp );
    TH1D *r_de2_wkp = new TH1D( "r_de2_wkp", "", bins_momega, h_momega_ll, h_momega_ul );
    norm_res( h_de_mw_proj2_wkp, p_de_mw_proj2_wkp, r_de2_wkp );
    TH1D *r_de3_wkp = new TH1D( "r_de3_wkp", "", bins_momega, h_momega_ll, h_momega_ul );
    norm_res( h_de_mw_proj3_wkp, p_de_mw_proj3_wkp, r_de3_wkp );
    TH1D *r_de4_wkp = new TH1D( "r_de4_wkp", "", bins_momega, h_momega_ll, h_momega_ul );
    norm_res( h_de_mw_proj4_wkp, p_de_mw_proj4_wkp, r_de4_wkp );
    TH1D *r_de5_wkp = new TH1D( "r_de5_wkp", "", bins_momega, h_momega_ll, h_momega_ul );
    norm_res( h_de_mw_proj5_wkp, p_de_mw_proj5_wkp, r_de5_wkp );    
    
  // Plot Mbc
  TH1D *h_mbc_svd1 = new TH1D( "h_mbc_svd1", "", bins_mbc, h_mbc_ll, h_mbc_ul);
  TH1D *h_mbc_svd2 = new TH1D( "h_mbc_svd2", "", bins_mbc, h_mbc_ll, h_mbc_ul);
  for( unsigned int i=0; i<data10.size(); i++ ){
    if( data10[i].GetExpNo() < 29 )
      h_mbc_svd1->Fill( data10[i].GetMbc());
    else
      h_mbc_svd2->Fill( data10[i].GetMbc());
  }
   
  TH1D *h_mbc_svd1_pdf = new TH1D( "h_mbc_svd1_pdf", "", bins_mbc, h_mbc_ll, h_mbc_ul );    
  TH1D *h_mbc_svd2_pdf = new TH1D( "h_mbc_svd2_pdf", "", bins_mbc, h_mbc_ll, h_mbc_ul );
   
   const double mbc_width =
    ( h_mbc_ul - h_mbc_ll ) / static_cast<double>(bins_mbc);
    
   for( unsigned int k=0; k < bins_mbc; k++ )
      {
     	 const double h_mbc_svd1 =
		h_mbc_ll + ((0.5+static_cast<double>(k))*mbc_width);     	 
	 const double h_mbc_svd2 =
		h_mbc_ll + ((0.5+static_cast<double>(k))*mbc_width);
     	 double h_pdf_mbc_svd1 = 0.0;
     	 double h_pdf_mbc_svd2 = 0.0;
     
	for( unsigned int l=0; l<data10.size(); l++ ){	
	      if( data10[l].GetExpNo() < 29 )
		h_pdf_mbc_svd1 += fcn.wksMBC(data10[l], h_mbc_svd1, par)/fcn.NormwksMBC("svd1", par);
	      else
 		h_pdf_mbc_svd2 += fcn.wksMBC(data10[l], h_mbc_svd2, par)/fcn.NormwksMBC("svd2", par);
	 }
      h_mbc_svd1_pdf->Fill(h_mbc_svd1, h_pdf_mbc_svd1);
      h_mbc_svd2_pdf->Fill(h_mbc_svd2, h_pdf_mbc_svd2);
      }
 
  h_mbc_svd1_pdf->Scale(mbc_width);
  h_mbc_svd2_pdf->Scale(mbc_width);  
  
  TH1D *h_mbc_svd1_wkp = new TH1D( "h_mbc_svd1_wkp", "", bins_mbc, h_mbc_ll, h_mbc_ul);
  TH1D *h_mbc_svd2_wkp = new TH1D( "h_mbc_svd2_wkp", "", bins_mbc, h_mbc_ll, h_mbc_ul);
  for( unsigned int i=0; i<data11.size(); i++ ){
    if( data11[i].GetExpNo() < 29 )
      h_mbc_svd1_wkp->Fill( data11[i].GetMbc());
    else
      h_mbc_svd2_wkp->Fill( data11[i].GetMbc());
  }
  
  TH1D *h_mbc_svd1_pdf_wkp = new TH1D( "h_mbc_svd1_pdf_wkp", "", bins_mbc, h_mbc_ll, h_mbc_ul );    
  TH1D *h_mbc_svd2_pdf_wkp = new TH1D( "h_mbc_svd2_pdf_wkp", "", bins_mbc, h_mbc_ll, h_mbc_ul );
    
   for( unsigned int k=0; k < bins_mbc; k++ )
      {
     	 const double h_mbc_svd1_wkp =
		h_mbc_ll + ((0.5+static_cast<double>(k))*mbc_width);     	 
	 const double h_mbc_svd2_wkp =
		h_mbc_ll + ((0.5+static_cast<double>(k))*mbc_width);
     	 double h_pdf_mbc_svd1_wkp = 0.0;
     	 double h_pdf_mbc_svd2_wkp = 0.0;
     
	for( unsigned int l=0; l<data11.size(); l++ ){	
	      if( data11[l].GetExpNo() < 29 )
		h_pdf_mbc_svd1_wkp += fcn.wkpMBC(data11[l], h_mbc_svd1_wkp, par)/fcn.NormwkpMBC("svd1", par);
	      else
 		h_pdf_mbc_svd2_wkp += fcn.wkpMBC(data11[l], h_mbc_svd2_wkp, par)/fcn.NormwkpMBC("svd2", par);
	 }
      h_mbc_svd1_pdf_wkp->Fill(h_mbc_svd1_wkp, h_pdf_mbc_svd1_wkp);
      h_mbc_svd2_pdf_wkp->Fill(h_mbc_svd2_wkp, h_pdf_mbc_svd2_wkp);
      }
 
  h_mbc_svd1_pdf_wkp->Scale(mbc_width);
  h_mbc_svd2_pdf_wkp->Scale(mbc_width);  
  
  // Plot omega h
  TH1D *h_oh_svd1 = new TH1D( "h_oh_svd1", "", bins_oh, h_homega_ll, h_homega_ul);
  TH1D *h_oh_svd2 = new TH1D( "h_oh_svd2", "", bins_oh, h_homega_ll, h_homega_ul);
  for( unsigned int i=0; i<data10.size(); i++ ){
    if( data10[i].GetExpNo() < 29 )
      h_oh_svd1->Fill( data10[i].GetHomega());
    else
      h_oh_svd2->Fill( data10[i].GetHomega());
  }
 
  
  TH1D *h_oh_svd1_pdf = new TH1D( "h_oh_svd1_pdf", "", bins_oh, h_homega_ll, h_homega_ul );    
  TH1D *h_oh_svd2_pdf = new TH1D( "h_oh_svd2_pdf", "", bins_oh, h_homega_ll, h_homega_ul );
   
   const double oh_width =
    ( h_homega_ul - h_homega_ll ) / static_cast<double>(bins_oh);
    
   for( unsigned int k=0; k < bins_oh; k++ )
      {
     	 const double h_oh_svd1 =
		h_homega_ll + ((0.5+static_cast<double>(k))*oh_width);     	 
	 const double h_oh_svd2 =
		h_homega_ll + ((0.5+static_cast<double>(k))*oh_width);
     	 double h_pdf_oh_svd1 = 0.0;
     	 double h_pdf_oh_svd2 = 0.0;
     
	for( unsigned int l=0; l<data10.size(); l++ ){	
	      if( data10[l].GetExpNo() < 29 )
		h_pdf_oh_svd1 += fcn.wksHw(data10[l], h_oh_svd1, par)/fcn.NormwksHw("svd1", par);
	      else
 		h_pdf_oh_svd2 += fcn.wksHw(data10[l], h_oh_svd2, par)/fcn.NormwksHw("svd2", par);
	 }
      h_oh_svd1_pdf->Fill(h_oh_svd1, h_pdf_oh_svd1);
      h_oh_svd2_pdf->Fill(h_oh_svd2, h_pdf_oh_svd2);
      }
 
  h_oh_svd1_pdf->Scale(oh_width);
  h_oh_svd2_pdf->Scale(oh_width);  
  
  TH1D *h_oh_svd1_wkp = new TH1D( "h_oh_svd1_wkp", "", bins_oh, h_homega_ll, h_homega_ul);
  TH1D *h_oh_svd2_wkp = new TH1D( "h_oh_svd2_wkp", "", bins_oh, h_homega_ll, h_homega_ul);
  for( unsigned int i=0; i<data11.size(); i++ ){
    if( data11[i].GetExpNo() < 29 )
      h_oh_svd1_wkp->Fill( data11[i].GetHomega());
    else
      h_oh_svd2_wkp->Fill( data11[i].GetHomega());
  }
  
  TH1D *h_oh_svd1_pdf_wkp = new TH1D( "h_oh_svd1_pdf_wkp", "", bins_oh, h_homega_ll, h_homega_ul );    
  TH1D *h_oh_svd2_pdf_wkp = new TH1D( "h_oh_svd2_pdf_wkp", "", bins_oh, h_homega_ll, h_homega_ul );
    
   for( unsigned int k=0; k < bins_oh; k++ )
      {
     	 const double h_oh_svd1_wkp =
		h_homega_ll + ((0.5+static_cast<double>(k))*oh_width);     	 
	 const double h_oh_svd2_wkp =
		h_homega_ll + ((0.5+static_cast<double>(k))*oh_width);
     	 double h_pdf_oh_svd1_wkp = 0.0;
     	 double h_pdf_oh_svd2_wkp = 0.0;
     
	for( unsigned int l=0; l<data11.size(); l++ ){	
	      if( data11[l].GetExpNo() < 29 )
		h_pdf_oh_svd1_wkp += fcn.wksHw(data11[l], h_oh_svd1_wkp, par)/fcn.NormwksHw("svd1", par);
	      else
 		h_pdf_oh_svd2_wkp += fcn.wksHw(data11[l], h_oh_svd2_wkp, par)/fcn.NormwksHw("svd2", par);
	 }
      h_oh_svd1_pdf_wkp->Fill(h_oh_svd1_wkp, h_pdf_oh_svd1_wkp);
      h_oh_svd2_pdf_wkp->Fill(h_oh_svd2_wkp, h_pdf_oh_svd2_wkp);
      }
 
  h_oh_svd1_pdf_wkp->Scale(oh_width);
  h_oh_svd2_pdf_wkp->Scale(oh_width);
 
  // Plot dt
  TH1D *h_dtb0_pdf_svd1 =
    new TH1D( "h_dtb0_pdf_svd1", "", bins_dt, h_dt_ll, h_dt_ul );  
  TH1D *h_dtb0_pdf_svd2 =
    new TH1D( "h_dtb0_pdf_svd2", "", bins_dt, h_dt_ll, h_dt_ul ); 
  TH1D *h_dtb0b_pdf_svd1 =
    new TH1D( "h_dtb0b_pdf_svd1", "", bins_dt, h_dt_ll, h_dt_ul );  
  TH1D *h_dtb0b_pdf_svd2 =
    new TH1D( "h_dtb0b_pdf_svd2", "", bins_dt, h_dt_ll, h_dt_ul );
    
  for( unsigned int l=0; l < bins_dt; l++ )
    {
	const double dt_svd1 = h_dtb0_pdf_svd1->GetBinCenter(l+1);
	const double dt_svd2 = h_dtb0_pdf_svd2->GetBinCenter(l+1);
     	double h_b0_pdf_svd1 = 0.0, h_b0_pdf_svd2 = 0.0, h_b0b_pdf_svd1 = 0.0, h_b0b_pdf_svd2 = 0.0;

     	 //Sum over all events
      	for( unsigned int m=0; m<data10.size(); m++ ){
	    if( data10[m].GetExpNo() < 29 ){
		h_b0_pdf_svd1 += fcn.wksDtQ( data10[m], dt_svd1, +1, par);
		h_b0b_pdf_svd1 += fcn.wksDtQ( data10[m], dt_svd1, -1, par );
	    }
	    else{
		h_b0_pdf_svd2 += fcn.wksDtQ( data10[m], dt_svd2, +1, par);
		h_b0b_pdf_svd2 += fcn.wksDtQ( data10[m], dt_svd2, -1, par);	    
	    }
	}
     	 h_dtb0_pdf_svd1->Fill(dt_svd1,h_b0_pdf_svd1);
     	 h_dtb0b_pdf_svd1->Fill(dt_svd1,h_b0b_pdf_svd1);     	
	 h_dtb0_pdf_svd2->Fill(dt_svd2,h_b0_pdf_svd2);
     	 h_dtb0b_pdf_svd2->Fill(dt_svd2,h_b0b_pdf_svd2);
    }
	 h_dtb0_pdf_svd1->Scale(h_dtb0_pdf_svd1->GetBinWidth(0));
	 h_dtb0b_pdf_svd1->Scale(h_dtb0b_pdf_svd1->GetBinWidth(0));	 
	 h_dtb0_pdf_svd2->Scale(h_dtb0_pdf_svd2->GetBinWidth(0));
	 h_dtb0b_pdf_svd2->Scale(h_dtb0b_pdf_svd2->GetBinWidth(0));

  TH1D *h_dt_b0_svd1 = new TH1D( "h_dt_b0_svd1", "", bins_dt, h_dt_ll, h_dt_ul);
  TH1D *h_dt_b0_svd2 = new TH1D( "h_dt_b0_svd2", "", bins_dt, h_dt_ll, h_dt_ul);
  TH1D *h_dt_b0b_svd1 = new TH1D( "h_dt_b0b_svd1", "", bins_dt, h_dt_ll, h_dt_ul);
  TH1D *h_dt_b0b_svd2 = new TH1D( "h_dt_b0b_svd2", "", bins_dt, h_dt_ll, h_dt_ul);

  for( unsigned int i=0; i<data10.size(); i++ ){
    if( data10[i].GetQ() == 1 ){
      if( data10[i].GetExpNo() < 29 )
	 h_dt_b0_svd1->Fill( data10[i].GetDeltat());	
      else 
	 h_dt_b0_svd2->Fill( data10[i].GetDeltat());		
    }
    else{	
      if( data10[i].GetExpNo() < 29 )
	 h_dt_b0b_svd1->Fill( data10[i].GetDeltat());	
      else 
	 h_dt_b0b_svd2->Fill( data10[i].GetDeltat());		
    }  
  }

  TH1D *h_dtb0_pdf_svd1_wkp =
    new TH1D( "h_dtb0_pdf_svd1_wkp ", "", bins_dt, h_dt_ll, h_dt_ul );  
  TH1D *h_dtb0_pdf_svd2_wkp  =
    new TH1D( "h_dtb0_pdf_svd2_wkp ", "", bins_dt, h_dt_ll, h_dt_ul ); 
  TH1D *h_dtb0b_pdf_svd1_wkp  =
    new TH1D( "h_dtb0b_pdf_svd1_wkp ", "", bins_dt, h_dt_ll, h_dt_ul );  
  TH1D *h_dtb0b_pdf_svd2_wkp  =
    new TH1D( "h_dtb0b_pdf_svd2_wkp ", "", bins_dt, h_dt_ll, h_dt_ul );
    
  for( unsigned int l=0; l < bins_dt; l++ )
    {
	const double dt_svd1 = h_dtb0_pdf_svd1_wkp ->GetBinCenter(l+1);
	const double dt_svd2 = h_dtb0_pdf_svd2_wkp ->GetBinCenter(l+1);
     	double h_b0_pdf_svd1 = 0.0, h_b0_pdf_svd2 = 0.0, h_b0b_pdf_svd1 = 0.0, h_b0b_pdf_svd2 = 0.0;

     	 //Sum over all events
      	for( unsigned int m=0; m<data11.size(); m++ ){
	    if( data11[m].GetExpNo() < 29 ){
		h_b0_pdf_svd1 += fcn.wkpDtQ( data11[m], dt_svd1, +1, par);
		h_b0b_pdf_svd1 += fcn.wkpDtQ( data11[m], dt_svd1, -1, par );
	    }
	    else{
		h_b0_pdf_svd2 += fcn.wkpDtQ( data11[m], dt_svd2, +1, par);
		h_b0b_pdf_svd2 += fcn.wkpDtQ( data11[m], dt_svd2, -1, par);	    
	    }
	}
     	 h_dtb0_pdf_svd1_wkp ->Fill(dt_svd1,h_b0_pdf_svd1);
     	 h_dtb0b_pdf_svd1_wkp ->Fill(dt_svd1,h_b0b_pdf_svd1);     	
	 h_dtb0_pdf_svd2_wkp ->Fill(dt_svd2,h_b0_pdf_svd2);
     	 h_dtb0b_pdf_svd2_wkp ->Fill(dt_svd2,h_b0b_pdf_svd2);
    }
	 h_dtb0_pdf_svd1_wkp ->Scale(h_dtb0_pdf_svd1->GetBinWidth(0));
	 h_dtb0b_pdf_svd1_wkp ->Scale(h_dtb0b_pdf_svd1->GetBinWidth(0));	 
	 h_dtb0_pdf_svd2_wkp ->Scale(h_dtb0_pdf_svd2->GetBinWidth(0));
	 h_dtb0b_pdf_svd2_wkp ->Scale(h_dtb0b_pdf_svd2->GetBinWidth(0));

  TH1D *h_dt_b0_svd1_wkp  = new TH1D( "h_dt_b0_svd1_wkp ", "", bins_dt, h_dt_ll, h_dt_ul);
  TH1D *h_dt_b0_svd2_wkp  = new TH1D( "h_dt_b0_svd2_wkp ", "", bins_dt, h_dt_ll, h_dt_ul);
  TH1D *h_dt_b0b_svd1_wkp  = new TH1D( "h_dt_b0b_svd1_wkp ", "", bins_dt, h_dt_ll, h_dt_ul);
  TH1D *h_dt_b0b_svd2_wkp  = new TH1D( "h_dt_b0b_svd2_wkp ", "", bins_dt, h_dt_ll, h_dt_ul);

  for( unsigned int i=0; i<data11.size(); i++ ){
    if( data11[i].GetQ() == 1 ){
      if( data11[i].GetExpNo() < 29 )
	 h_dt_b0_svd1_wkp ->Fill( data11[i].GetDeltat());	
      else 
	 h_dt_b0_svd2_wkp ->Fill( data11[i].GetDeltat());		
    }
    else{	
      if( data11[i].GetExpNo() < 29 )
	 h_dt_b0b_svd1_wkp ->Fill( data11[i].GetDeltat());	
      else 
	 h_dt_b0b_svd2_wkp ->Fill( data11[i].GetDeltat());		
    }  
  }

  // Plot fd
  TH1D *h_fd_svd1 = new TH1D( "h_fd_svd1", "", 2*bins_fd, h_fd_ll, h_fd_ul);
  TH1D *h_fd_svd2 = new TH1D( "h_fd_svd2", "", 2*bins_fd, h_fd_ll, h_fd_ul);
  for( unsigned int i=0; i<data10.size(); i++ ){
//if( 0.875 < data10[i].GetR()&& data10[i].GetR() <= 1.){
   if( data10[i].GetExpNo() < 29 )
    h_fd_svd1->Fill( data10[i].GetFD());
   else
    h_fd_svd2->Fill( data10[i].GetFD());
//}
  }
  TH1D *h_fd_pdf_svd1 =
    new TH1D( "h_fd_pdf_svd1", "", 2*bins_fd, h_fd_ll, h_fd_ul );  
  TH1D *h_fd_pdf_svd2 =
    new TH1D( "h_fd_pdf_svd2", "", 2*bins_fd, h_fd_ll, h_fd_ul );
    
   const double fd_width =
    ( h_fd_ul - h_fd_ll ) / static_cast<double>(2*bins_fd);
    
   for( unsigned int k=0; k < 2*bins_fd; k++ )
      {
     	 const double h_fd_svd1 =
		h_fd_ll + ((0.5+static_cast<double>(k))*fd_width);     	 
	 const double h_fd_svd2 =
		h_fd_ll + ((0.5+static_cast<double>(k))*fd_width);
     	 double h_pdf_fd_svd1 = 0.0;
     	 double h_pdf_fd_svd2 = 0.0;
     
	for( unsigned int l=0; l<data10.size(); l++ ){
 	  int r_bin = set_rbin( data10[l].GetR() );
	      if( data10[l].GetExpNo() < 29 )
		h_pdf_fd_svd1 += fcn.wksFD(data10[l], h_fd_svd1, r_bin, par)/fcn.NormwksFD("svd1", r_bin, par);
		else
 		h_pdf_fd_svd2 += fcn.wksFD(data10[l], h_fd_svd2, r_bin, par)/fcn.NormwksFD("svd2", r_bin, par);
	  }
	h_fd_pdf_svd1->Fill(h_fd_svd1, h_pdf_fd_svd1);
	h_fd_pdf_svd2->Fill(h_fd_svd2, h_pdf_fd_svd2);
      }
 
  h_fd_pdf_svd1->Scale(fd_width);
  h_fd_pdf_svd2->Scale(fd_width);  

  TH1D *h_fd_svd1_wkp = new TH1D( "h_fd_svd1_wkp", "", 2*bins_fd, h_fd_ll, h_fd_ul);
  TH1D *h_fd_svd2_wkp = new TH1D( "h_fd_svd2_wkp", "", 2*bins_fd, h_fd_ll, h_fd_ul);
  for( unsigned int i=0; i<data11.size(); i++ ){
   if( data11[i].GetExpNo() < 29 )
    h_fd_svd1_wkp->Fill( data11[i].GetFD());
   else
    h_fd_svd2_wkp->Fill( data11[i].GetFD());
  }
  TH1D *h_fd_pdf_svd1_wkp =
    new TH1D( "h_fd_pdf_svd1_wkp", "", 2*bins_fd, h_fd_ll, h_fd_ul );  
  TH1D *h_fd_pdf_svd2_wkp=
    new TH1D( "h_fd_pdf_svd2_wkp", "", 2*bins_fd, h_fd_ll, h_fd_ul );
    
   for( unsigned int k=0; k < 2*bins_fd; k++ )
      {
     	 const double h_fd_svd1_wkp =
		h_fd_ll + ((0.5+static_cast<double>(k))*fd_width);     	 
	 const double h_fd_svd2_wkp =
		h_fd_ll + ((0.5+static_cast<double>(k))*fd_width);
     	 double h_pdf_fd_svd1_wkp = 0.0;
     	 double h_pdf_fd_svd2_wkp = 0.0;
     
	for( unsigned int l=0; l<data11.size(); l++ ){
 	  int r_bin = set_rbin( data11[l].GetR() );
	      if( data11[l].GetExpNo() < 29 )
		h_pdf_fd_svd1_wkp += fcn.wkpFD(data11[l], h_fd_svd1_wkp, r_bin, par)/fcn.NormwkpFD("svd1", r_bin, par);
		else
 		h_pdf_fd_svd2_wkp += fcn.wkpFD(data11[l], h_fd_svd2_wkp, r_bin, par)/fcn.NormwkpFD("svd2", r_bin, par);
	  }
	
	h_fd_pdf_svd1_wkp->Fill(h_fd_svd1_wkp, h_pdf_fd_svd1_wkp);
	h_fd_pdf_svd2_wkp->Fill(h_fd_svd2_wkp, h_pdf_fd_svd2_wkp);
      }
 
  h_fd_pdf_svd1_wkp->Scale(fd_width);
  h_fd_pdf_svd2_wkp->Scale(fd_width); 
  
  TH1D *r_de = new TH1D( "r_de", "", bins_de, h_de_ll, h_de_ul );
    norm_res( h_de_sig, h_de_pdf, r_de );   
  TH1D *r_de_wkp = new TH1D( "r_de_wkp", "", bins_de, h_de_ll, h_de_ul );
    norm_res( h_de_sig_wkp, h_de_pdf_wkp, r_de_wkp );    
  TH1D *r_om = new TH1D( "r_om", "", bins_momega, h_momega_ll, h_momega_ul );
    norm_res( h_om_sig, h_om_pdf, r_om );   
  TH1D *r_om_wkp = new TH1D( "r_om_wkp", "", bins_momega, h_momega_ll, h_momega_ul );
    norm_res( h_om_sig_wkp, h_om_pdf_wkp, r_om_wkp );  
 /* TH1D *r_de_svd2 = new TH1D( "r_de_svd2", "", bins_de, h_de_ll, h_de_ul );
    norm_res( h_de_svd2, h_de_svd2_pdf, r_de_svd2 );  
  TH1D *r_de_svd1_wkp = new TH1D( "r_de_svd1_wkp", "", bins_de, h_de_ll, h_de_ul );
    norm_res( h_de_svd1_wkp, h_de_svd1_pdf_wkp, r_de_svd1_wkp);  
  TH1D *r_de_svd2_wkp = new TH1D( "r_de_svd2_wkp", "", bins_de, h_de_ll, h_de_ul );
    norm_res( h_de_svd2_wkp, h_de_svd2_pdf_wkp, r_de_svd2_wkp );  */
  TH1D *r_oh1 = new TH1D( "r_oh1", "", bins_oh, h_homega_ll, h_homega_ul );
    norm_res( h_oh_svd1, h_oh_svd1_pdf, r_oh1 );  
  TH1D *r_oh2 = new TH1D( "r_oh2", "", bins_oh, h_homega_ll, h_homega_ul );
    norm_res( h_oh_svd2, h_oh_svd2_pdf, r_oh2 );  
  TH1D *r_oh1_wkp = new TH1D( "r_oh1_wkp", "", bins_oh, h_homega_ll, h_homega_ul );
    norm_res( h_oh_svd1_wkp, h_oh_svd1_pdf_wkp, r_oh1_wkp);  
  TH1D *r_oh2_wkp = new TH1D( "r_oh2_wkp", "", bins_oh, h_homega_ll, h_homega_ul );
    norm_res( h_oh_svd2_wkp, h_oh_svd2_pdf_wkp, r_oh2_wkp );  
/*  TH1D *r_om1 = new TH1D( "r_om1", "", bins_om, h_momega_ll, h_momega_ul );
    norm_res( h_om_svd1, h_om_svd1_pdf, r_om1 );  
  TH1D *r_om2 = new TH1D( "r_om2", "", bins_om, h_momega_ll, h_momega_ul );
    norm_res( h_om_svd2, h_om_svd2_pdf, r_om2 );  
   TH1D *r_om1_wkp = new TH1D( "r_om1_wkp", "", bins_om, h_momega_ll, h_momega_ul );
    norm_res( h_om_svd1_wkp, h_om_svd1_pdf_wkp, r_om1_wkp);  
  TH1D *r_om2_wkp = new TH1D( "r_om2_wkp", "", bins_om, h_momega_ll, h_momega_ul );
    norm_res( h_om_svd2_wkp, h_om_svd2_pdf_wkp, r_om2_wkp );  */
  TH1D *r_fd1 = new TH1D( "r_fd1", "", 2*bins_fd, h_fd_ll, h_fd_ul );
    norm_res( h_fd_svd1, h_fd_pdf_svd1, r_fd1 );  
  TH1D *r_fd2 = new TH1D( "r_fd2", "", 2*bins_fd, h_fd_ll, h_fd_ul );
    norm_res( h_fd_svd2, h_fd_pdf_svd2, r_fd2 );  
  TH1D *r_fd1_wkp = new TH1D( "r_fd1_wkp", "", 2*bins_fd, h_fd_ll, h_fd_ul );
    norm_res( h_fd_svd1_wkp, h_fd_pdf_svd1_wkp, r_fd1_wkp);  
  TH1D *r_fd2_wkp = new TH1D( "r_fd2_wkp", "", 2*bins_fd, h_fd_ll, h_fd_ul );
    norm_res( h_fd_svd2_wkp, h_fd_pdf_svd2_wkp, r_fd2_wkp );    
  TH1D *r_mbc1 = new TH1D( "r_mbc1", "", bins_mbc, h_mbc_ll, h_mbc_ul );
    norm_res( h_mbc_svd1, h_mbc_svd1_pdf, r_mbc1 );  
  TH1D *r_mbc2 = new TH1D( "r_mbc2", "", bins_mbc, h_mbc_ll, h_mbc_ul );
    norm_res( h_mbc_svd2, h_mbc_svd2_pdf, r_mbc2 );  
  TH1D *r_mbc1_wkp = new TH1D( "r_mbc1_wkp", "", bins_mbc, h_mbc_ll, h_mbc_ul );
    norm_res( h_mbc_svd1_wkp, h_mbc_svd1_pdf_wkp, r_mbc1_wkp);  
  TH1D *r_mbc2_wkp = new TH1D( "r_mbc2_wkp", "", bins_mbc, h_mbc_ll, h_mbc_ul );
    norm_res( h_mbc_svd2_wkp, h_mbc_svd2_pdf_wkp, r_mbc2_wkp );     
/*    
  const double chisq_de_svd1 = get_chisq(h_de_svd1, h_de_svd1_pdf, 8);
  std::stringstream ss_chisq_de_svd1;
  ss_chisq_de_svd1 << std::setprecision(2) << chisq_de_svd1;
  const std::string schisq_de_svd1 = "#scale[1.5]{#chi^{2}/ndf=" + ss_chisq_de_svd1.str() + "}";
  
  const double chisq_de_svd2 = get_chisq(h_de_svd2, h_de_svd2_pdf, 8);
  std::stringstream ss_chisq_de_svd2;
  ss_chisq_de_svd2 << std::setprecision(2) << chisq_de_svd2;
  const std::string schisq_de_svd2 = "#scale[1.5]{#chi^{2}/ndf=" + ss_chisq_de_svd2.str() + "}";  

  TCanvas *canvas11 = new TCanvas("canvas11","",1024, 768);
  normres_style(canvas11);
  canvas11->cd(1);
  h_de_svd1->Draw("E1");
  h_de_svd1_pdf->SetLineColor(kBlue);
  h_de_svd1_pdf->Draw("C SAME");
  
 TLatex *ltxt_de_svd1 = new TLatex(h_de_ll+0.005, 0.9*h_de_svd1->GetMaximum(),
			    schisq_de_svd1.c_str());
  ltxt_de_svd1->Draw();
  
  canvas11->cd(2);
  normres_style(r_de1, normres_sigma, "B^{0} #DeltaE", "GeV");
  r_de1->Draw("E1");
  for( unsigned int i=0; i<normres_sigma.size(); i++ )
    normres_sigma.at(i).Draw();  
  canvas11->Update();
  canvas11->Print("de_svd1_wks.eps");
  
  TCanvas *canvas12 = new TCanvas("canvas12","",1024, 768);
  normres_style(canvas12);
  canvas12->cd(1);
  h_de_svd2->Draw("E1");
  h_de_svd2_pdf->SetLineColor(kBlue);
  h_de_svd2_pdf->Draw("C SAME");
  
  TLatex *ltxt_de_svd2 = new TLatex(h_de_ll+0.005, 0.9*h_de_svd2->GetMaximum(),
			    schisq_de_svd2.c_str());
  ltxt_de_svd2->Draw();
  
  canvas12->cd(2);
  normres_style(r_de2, normres_sigma, "B^{0} #DeltaE", "GeV");
  r_de2->Draw("E1");
  for( unsigned int i=0; i<normres_sigma.size(); i++ )
    normres_sigma.at(i).Draw();  
  canvas12->Update();
  canvas12->Print("de_svd2_wks.eps");
  
  const double chisq_de_svd1_wkp = get_chisq(h_de_svd1_wkp, h_de_svd1_pdf_wkp, 8);
  std::stringstream ss_chisq_de_svd1_wkp;
  ss_chisq_de_svd1_wkp << std::setprecision(2) << chisq_de_svd1_wkp;
  const std::string schisq_de_svd1_wkp = "#scale[1.5]{#chi^{2}/ndf=" + ss_chisq_de_svd1_wkp.str() + "}";
  
  const double chisq_de_svd2_wkp = get_chisq(h_de_svd2_wkp, h_de_svd2_pdf_wkp, 8);
  std::stringstream ss_chisq_de_svd2_wkp;
  ss_chisq_de_svd2_wkp << std::setprecision(2) << chisq_de_svd2_wkp;
  const std::string schisq_de_svd2_wkp = "#scale[1.5]{#chi^{2}/ndf=" + ss_chisq_de_svd2_wkp.str() + "}";  
  
  TCanvas *canvas11_wkp = new TCanvas("canvas11_wkp","",1024, 768);
  normres_style(canvas11_wkp);
  canvas11_wkp->cd(1);
  h_de_svd1_wkp->Draw("E1");
  h_de_svd1_pdf_wkp->SetLineColor(kBlue);
  h_de_svd1_pdf_wkp->Draw("C SAME");
  
 TLatex *ltxt_de_svd1_wkp = new TLatex(h_de_ll+0.005, 0.9*h_de_svd1_wkp->GetMaximum(),
			    schisq_de_svd1_wkp.c_str());
  ltxt_de_svd1_wkp->Draw();
  
  canvas11_wkp->cd(2);
  normres_style(r_de1, normres_sigma, "B^{+} #DeltaE", "GeV");
  r_de1_wkp->Draw("E1");
  for( unsigned int i=0; i<normres_sigma.size(); i++ )
    normres_sigma.at(i).Draw();  
  canvas11_wkp->Update();
  canvas11_wkp->Print("de_svd1_wkp.eps");
  
  TCanvas *canvas12_wkp = new TCanvas("canvas12_wkp","",1024, 768);
  normres_style(canvas12_wkp);
  canvas12_wkp->cd(1);
  h_de_svd2_wkp->Draw("E1");
  h_de_svd2_pdf_wkp->SetLineColor(kBlue);
  h_de_svd2_pdf_wkp->Draw("C SAME");
  
  TLatex *ltxt_de_svd2_wkp = new TLatex(h_de_ll+0.005, 0.9*h_de_svd2_wkp->GetMaximum(),
			    schisq_de_svd2_wkp.c_str());
  ltxt_de_svd2_wkp->Draw();
  
  canvas12_wkp->cd(2);
  normres_style(r_de2, normres_sigma, "B^{+} #DeltaE", "GeV");
  r_de2_wkp->Draw("E1");
  for( unsigned int i=0; i<normres_sigma.size(); i++ )
    normres_sigma.at(i).Draw();  
  canvas12_wkp->Update();
  canvas12_wkp->Print("de_svd2_wkp.eps");

  
  TCanvas *canvas101 = new TCanvas("canvas101","",1024, 1024);

  hist_style(h_de_svd1_wkp, "#DeltaE", "GeV");
  //h_de_svd1_wkp>SetMinimum(0.01);
  h_de_svd1_wkp->Draw("E1");
  h_de_svd1->Draw("E1 SAME");
  h_de_svd1_pdf->SetLineColor(kBlue);
  h_de_svd1_pdf->Draw("C SAME");  
  h_de_svd1_pdf_wkp->SetLineColor(kRed);
  h_de_svd1_pdf_wkp->Draw("C SAME");
  
  canvas101->Update();
  canvas101->Print("de_svd1.eps");
  
  TCanvas *canvas102 = new TCanvas("canvas102","",1024, 768);
  hist_style(h_de_svd2_wkp, "#DeltaE", "GeV");
  h_de_svd2_wkp->Draw("E1");
  h_de_svd2->Draw("E1");
  h_de_svd2_pdf->SetLineColor(kBlue);
  h_de_svd2_pdf->Draw("C SAME");  
  h_de_svd2_pdf_wkp->SetLineColor(kRed);
  h_de_svd2_pdf_wkp->Draw("C SAME");
  canvas102->Update();
  canvas102->Print("de_svd2.eps");
  

  const double chisq_om_svd1 = get_chisq(h_om_svd1, h_om_svd1_pdf, 4);
  std::stringstream ss_chisq_om_svd1;
  ss_chisq_om_svd1 << std::setprecision(2) << chisq_om_svd1;
  const std::string schisq_om_svd1 = "#scale[1.5]{#chi^{2}/ndf=" + ss_chisq_om_svd1.str() + "}";
  
  const double chisq_om_svd2 = get_chisq(h_om_svd2, h_om_svd2_pdf, 4);
  std::stringstream ss_chisq_om_svd2;
  ss_chisq_om_svd2 << std::setprecision(2) << chisq_om_svd2;
  const std::string schisq_om_svd2 = "#scale[1.5]{#chi^{2}/ndf=" + ss_chisq_om_svd2.str() + "}";  
  
  TCanvas *canvas21 = new TCanvas("canvas21","",1024, 768);
  normres_style(canvas21);
  canvas21->cd(1);
  h_om_svd1->Draw("E1");
  h_om_svd1_pdf->SetLineColor(kBlue);
  h_om_svd1_pdf->Draw("C SAME");
  TLatex *ltxt_om_svd1 = new TLatex(h_momega_ll+0.005, 0.9*h_om_svd1->GetMaximum(),
			    schisq_om_svd1.c_str());
  ltxt_om_svd1->Draw();  
  canvas21->cd(2);
  normres_style(r_om1, normres_sigma, "B^{0} H_{3#pi}", "");
  r_om1->Draw("E1");
  for( unsigned int i=0; i<normres_sigma.size(); i++ )
    normres_sigma.at(i).Draw();  
  canvas21->Update();
  canvas21->Print("omega_m_svd1_wks.eps");
  
  TCanvas *canvas22 = new TCanvas("canvas22","",1024, 768);
  normres_style(canvas22);
  canvas22->cd(1);
  h_om_svd2->Draw("E1");
  h_om_svd2_pdf->SetLineColor(kBlue);
  h_om_svd2_pdf->Draw("C SAME");
  TLatex *ltxt_om_svd2 = new TLatex(h_momega_ll+0.005, 0.9*h_om_svd2->GetMaximum(),
			    schisq_om_svd2.c_str());
  ltxt_om_svd2->Draw();
  canvas22->cd(2);
  normres_style(r_om1, normres_sigma, "B^{0} H_{3#pi}", "");
  r_om1->Draw("E1");
  for( unsigned int i=0; i<normres_sigma.size(); i++ )
    normres_sigma.at(i).Draw();  
  canvas22->Update();
  canvas22->Print("omega_m_svd2_wks.eps");

  const double chisq_om_svd1_wkp = get_chisq(h_om_svd1_wkp, h_om_svd1_pdf_wkp, 4);
  std::stringstream ss_chisq_om_svd1_wkp;
  ss_chisq_om_svd1_wkp << std::setprecision(2) << chisq_om_svd1_wkp;
  const std::string schisq_om_svd1_wkp = "#scale[1.5]{#chi^{2}/ndf=" + ss_chisq_om_svd1_wkp.str() + "}";
  
  const double chisq_om_svd2_wkp = get_chisq(h_om_svd2_wkp, h_om_svd2_pdf_wkp, 4);
  std::stringstream ss_chisq_om_svd2_wkp;
  ss_chisq_om_svd2_wkp << std::setprecision(2) << chisq_om_svd2_wkp;
  const std::string schisq_om_svd2_wkp= "#scale[1.5]{#chi^{2}/ndf=" + ss_chisq_om_svd2_wkp.str() + "}";  
  
  TCanvas *canvas21_wkp = new TCanvas("canvas21_wkp","",1024, 768);
  normres_style(canvas21_wkp);
  canvas21_wkp->cd(1);
  h_om_svd1_wkp->Draw("E1");
  h_om_svd1_pdf_wkp->SetLineColor(kBlue);
  h_om_svd1_pdf_wkp->Draw("C SAME");
  TLatex *ltxt_om_svd1_wkp = new TLatex(h_momega_ll+0.005, 0.9*h_om_svd1_wkp->GetMaximum(),
			    schisq_om_svd1_wkp.c_str());
  ltxt_om_svd1_wkp->Draw();  
  canvas21_wkp->cd(2);
  normres_style(r_om1_wkp, normres_sigma, "B^{+} H_{3#pi}", "");
  r_om1_wkp->Draw("E1");
  for( unsigned int i=0; i<normres_sigma.size(); i++ )
    normres_sigma.at(i).Draw();  
  canvas21_wkp->Update();
  canvas21_wkp->Print("omega_m_svd1_wkp.eps");
  
  TCanvas *canvas22_wkp = new TCanvas("canvas22_wkp","",1024, 768);
  normres_style(canvas22_wkp);
  canvas22_wkp->cd(1);
  h_om_svd2_wkp->Draw("E1");
  h_om_svd2_pdf_wkp->SetLineColor(kBlue);
  h_om_svd2_pdf_wkp->Draw("C SAME");
  TLatex *ltxt_om_svd2_wkp = new TLatex(h_momega_ll+0.005, 0.9*h_om_svd2_wkp->GetMaximum(),
			    schisq_om_svd2_wkp.c_str());
  ltxt_om_svd2_wkp->Draw();  
  canvas22_wkp->cd(2);
  normres_style(r_om2_wkp, normres_sigma, "B^{+} H_{3#pi}", "");
  r_om2_wkp->Draw("E1");
  for( unsigned int i=0; i<normres_sigma.size(); i++ )
    normres_sigma.at(i).Draw();  
  canvas22_wkp->Update();
  canvas22_wkp->Print("omega_m_svd2_wkp.eps");
 */ 

   
  const double chisq_de= get_chisq(h_de_sig, h_de_pdf, 10);
  std::stringstream ss_chisq_de;
  ss_chisq_de<< std::setprecision(2) << chisq_de;
  const std::string schisq_de = "#scale[1.5]{#chi^{2}/ndf=" + ss_chisq_de.str() + "}";  
  
  TCanvas *canvas1 = new TCanvas("canvas1","",1024, 768);
  normres_style(canvas1);
  canvas1->cd(1);
  hist_style(h_de_sig, "#Delta E", "GeV");
  h_de_sig->SetMinimum(0.01);
  h_de_sig->Draw("E1");
  h_de_pdf->SetLineColor(kBlue);
  h_de_pdf->Draw("HIST SAME");
    TLatex *ltxt_de = new TLatex(h_de_ll+0.005, 0.9*h_de_sig->GetMaximum(),
			    schisq_de.c_str());
  ltxt_de->Draw();
  canvas1->cd(2);
  normres_style(r_de, normres_sigma, "#Delta E", "GeV");
  r_de->Draw("E1");
  for( unsigned int i=0; i<normres_sigma.size(); i++ )
    normres_sigma.at(i).Draw();
  canvas1->Update();
  canvas1->Print("de_wks.eps");

  const double chisq_om= get_chisq(h_om_sig, h_om_pdf, 8);
  std::stringstream ss_chisq_om;
  ss_chisq_om<< std::setprecision(2) << chisq_om;
  const std::string schisq_om = "#scale[1.5]{#chi^{2}/ndf=" + ss_chisq_om.str() + "}";  
  
  TCanvas *canvas2 = new TCanvas("canvas2","",1024, 768);
  normres_style(canvas2);
  canvas2->cd(1);
  hist_style(h_de_sig, "#M_{#pi^{+}#pi^{-}#pi^{0}}", "GeV/c^{2}");
  h_om_sig->SetMinimum(0.01);
  h_om_sig->Draw("E1");
  h_om_pdf->SetLineColor(kBlue);
  h_om_pdf->Draw("HIST SAME");
    TLatex *ltxt_om = new TLatex(h_momega_ll+0.005, 0.9*h_om_sig->GetMaximum(),
			    schisq_om.c_str());
  ltxt_om->Draw();
  canvas2->cd(2);
  normres_style(r_om, normres_sigma, "#M_{#pi^{+}#pi^{-}#pi^{0}}", "GeV/c^{2}");
  r_om->Draw("E1");
  for( unsigned int i=0; i<normres_sigma.size(); i++ )
    normres_sigma.at(i).Draw();
  canvas2->Update();
  canvas2->Print("omega_m_wks.eps");

  const double chisq_de_wkp= get_chisq(h_de_sig_wkp, h_de_pdf_wkp, 10);
  std::stringstream ss_chisq_de_wkp;
  ss_chisq_de_wkp<< std::setprecision(2) << chisq_de_wkp;
  const std::string schisq_de_wkp = "#scale[1.5]{#chi^{2}/ndf=" + ss_chisq_de_wkp.str() + "}";  
  
  TCanvas *canvas1_wkp = new TCanvas("canvas1_wkp","",1024, 768);
  normres_style(canvas1_wkp);
  canvas1_wkp->cd(1);
  hist_style(h_de_sig_wkp, "#Delta E", "GeV");
  h_de_sig_wkp->SetMinimum(0.01);
  h_de_sig_wkp->Draw("E1");
  h_de_pdf_wkp->SetLineColor(kBlue);
  h_de_pdf_wkp->Draw("HIST SAME");
    TLatex *ltxt_de_wkp = new TLatex(h_de_ll+0.005, 0.9*h_de_sig_wkp->GetMaximum(),
			    schisq_de_wkp.c_str());
  ltxt_de_wkp->Draw();
  canvas1_wkp->cd(2);
  normres_style(r_de_wkp, normres_sigma, "#Delta E", "GeV");
  r_de_wkp->Draw("E1");
  for( unsigned int i=0; i<normres_sigma.size(); i++ )
    normres_sigma.at(i).Draw();
  canvas1_wkp->Update();
  canvas1_wkp->Print("de_wkp.eps");

  const double chisq_om_wkp= get_chisq(h_om_sig_wkp, h_om_pdf_wkp, 8);
  std::stringstream ss_chisq_om_wkp;
  ss_chisq_om_wkp<< std::setprecision(2) << chisq_om_wkp;
  const std::string schisq_om_wkp = "#scale[1.5]{#chi^{2}/ndf=" + ss_chisq_om_wkp.str() + "}";  
  
  TCanvas *canvas2_wkp = new TCanvas("canvas2_wkp","",1024, 768);
  normres_style(canvas2_wkp);
  canvas2_wkp->cd(1);
  hist_style(h_de_sig_wkp, "#M_{#pi^{+}#pi^{-}#pi^{0}}", "GeV/c^{2}");
  h_om_sig_wkp->SetMinimum(0.01);
  h_om_sig_wkp->Draw("E1");
  h_om_pdf_wkp->SetLineColor(kBlue);
  h_om_pdf_wkp->Draw("HIST SAME");
  TLatex *ltxt_om_wkp = new TLatex(h_momega_ll+0.005, 0.9*h_om_sig_wkp->GetMaximum(),
			    schisq_om_wkp.c_str());
  ltxt_om_wkp->Draw();
  canvas2_wkp->cd(2);
  normres_style(r_om_wkp, normres_sigma, "#M_{#pi^{+}#pi^{-}#pi^{0}}", "GeV/c^{2}");
  r_om_wkp->Draw("E1");
  for( unsigned int i=0; i<normres_sigma.size(); i++ )
    normres_sigma.at(i).Draw();
  canvas2_wkp->Update();
  canvas2_wkp->Print("omega_m_wkp.eps"); 
  
  const double chisq_oh_svd1 = get_chisq(h_oh_svd1, h_oh_svd1_pdf, 4);
  std::stringstream ss_chisq_oh_svd1;
  ss_chisq_oh_svd1 << std::setprecision(2) << chisq_oh_svd1;
  const std::string schisq_oh_svd1 = "#scale[1.5]{#chi^{2}/ndf=" + ss_chisq_oh_svd1.str() + "}";
  
  const double chisq_oh_svd2 = get_chisq(h_oh_svd2, h_oh_svd2_pdf, 4);
  std::stringstream ss_chisq_oh_svd2;
  ss_chisq_oh_svd2 << std::setprecision(2) << chisq_oh_svd2;
  const std::string schisq_oh_svd2 = "#scale[1.5]{#chi^{2}/ndf=" + ss_chisq_oh_svd2.str() + "}";  
  
  TCanvas *canvas31 = new TCanvas("canvas31","",1024, 768);
  normres_style(canvas31);
  canvas31->cd(1);
  h_oh_svd1->Draw("E1");
  h_oh_svd1_pdf->SetLineColor(kBlue);
  h_oh_svd1_pdf->Draw("C SAME");
  TLatex *ltxt_oh_svd1 = new TLatex(h_homega_ll+0.005, 0.9*h_oh_svd1->GetMaximum(),
			    schisq_oh_svd1.c_str());
  ltxt_oh_svd1->Draw();  
  canvas31->cd(2);
  normres_style(r_oh1, normres_sigma, "B^{0} H_{3#pi}", "");
  r_oh1->Draw("E1");
  for( unsigned int i=0; i<normres_sigma.size(); i++ )
    normres_sigma.at(i).Draw();  
  canvas31->Update();
  canvas31->Print("omega_h_svd1_wks.eps");
  
  TCanvas *canvas32 = new TCanvas("canvas32","",1024, 768);
  normres_style(canvas32);
  canvas32->cd(1);
  h_oh_svd2->Draw("E1");
  h_oh_svd2_pdf->SetLineColor(kBlue);
  h_oh_svd2_pdf->Draw("C SAME");
  TLatex *ltxt_oh_svd2 = new TLatex(h_homega_ll+0.005, 0.9*h_oh_svd2->GetMaximum(),
			    schisq_oh_svd2.c_str());
  ltxt_oh_svd2->Draw();
  canvas32->cd(2);
  normres_style(r_oh1, normres_sigma, "B^{0} H_{3#pi}", "");
  r_oh1->Draw("E1");
  for( unsigned int i=0; i<normres_sigma.size(); i++ )
    normres_sigma.at(i).Draw();  
  canvas32->Update();
  canvas32->Print("omega_h_svd2_wks.eps");

  const double chisq_oh_svd1_wkp = get_chisq(h_oh_svd1_wkp, h_oh_svd1_pdf_wkp, 4);
  std::stringstream ss_chisq_oh_svd1_wkp;
  ss_chisq_oh_svd1_wkp << std::setprecision(2) << chisq_oh_svd1_wkp;
  const std::string schisq_oh_svd1_wkp = "#scale[1.5]{#chi^{2}/ndf=" + ss_chisq_oh_svd1_wkp.str() + "}";
  
  const double chisq_oh_svd2_wkp = get_chisq(h_oh_svd2_wkp, h_oh_svd2_pdf_wkp, 4);
  std::stringstream ss_chisq_oh_svd2_wkp;
  ss_chisq_oh_svd2_wkp << std::setprecision(2) << chisq_oh_svd2_wkp;
  const std::string schisq_oh_svd2_wkp= "#scale[1.5]{#chi^{2}/ndf=" + ss_chisq_oh_svd2_wkp.str() + "}";  
  
  TCanvas *canvas31_wkp = new TCanvas("canvas31_wkp","",1024, 768);
  normres_style(canvas31_wkp);
  canvas31_wkp->cd(1);
  h_oh_svd1_wkp->Draw("E1");
  h_oh_svd1_pdf_wkp->SetLineColor(kBlue);
  h_oh_svd1_pdf_wkp->Draw("C SAME");
  TLatex *ltxt_oh_svd1_wkp = new TLatex(h_homega_ll+0.005, 0.9*h_oh_svd1_wkp->GetMaximum(),
			    schisq_oh_svd1_wkp.c_str());
  ltxt_oh_svd1_wkp->Draw();  
  canvas31_wkp->cd(2);
  normres_style(r_oh1_wkp, normres_sigma, "B^{+} H_{3#pi}", "");
  r_oh1_wkp->Draw("E1");
  for( unsigned int i=0; i<normres_sigma.size(); i++ )
    normres_sigma.at(i).Draw();  
  canvas31_wkp->Update();
  canvas31_wkp->Print("omega_h_svd1_wkp.eps");
  
  TCanvas *canvas32_wkp = new TCanvas("canvas32_wkp","",1024, 768);
  normres_style(canvas32_wkp);
  canvas32_wkp->cd(1);
  h_oh_svd2_wkp->Draw("E1");
  h_oh_svd2_pdf_wkp->SetLineColor(kBlue);
  h_oh_svd2_pdf_wkp->Draw("C SAME");
  TLatex *ltxt_oh_svd2_wkp = new TLatex(h_homega_ll+0.005, 0.9*h_oh_svd2_wkp->GetMaximum(),
			    schisq_oh_svd2_wkp.c_str());
  ltxt_oh_svd2_wkp->Draw();  
  canvas32_wkp->cd(2);
  normres_style(r_oh2_wkp, normres_sigma, "B^{+} H_{3#pi}", "");
  r_oh2_wkp->Draw("E1");
  for( unsigned int i=0; i<normres_sigma.size(); i++ )
    normres_sigma.at(i).Draw();  
  canvas32_wkp->Update();
  canvas32_wkp->Print("omega_h_svd2_wkp.eps");
  
  
  const double chisq_dt_b0_svd1 = get_chisq(h_dt_b0_svd1, h_dtb0_pdf_svd1, 3);
  std::stringstream ss_chisq_dt_b0_svd1;
  ss_chisq_dt_b0_svd1 << std::setprecision(2) << chisq_dt_b0_svd1;
  const std::string schisq_dt_b0_svd1 = "#scale[1.5]{#chi^{2}/ndf=" + ss_chisq_dt_b0_svd1.str() + "}";
  
  const double chisq_dt_b0_svd2 = get_chisq(h_dt_b0_svd2, h_dtb0_pdf_svd2, 3);
  std::stringstream ss_chisq_dt_b0_svd2;
  ss_chisq_dt_b0_svd2 << std::setprecision(2) << chisq_dt_b0_svd2;
  const std::string schisq_dt_b0_svd2 = "#scale[1.5]{#chi^{2}/ndf=" + ss_chisq_dt_b0_svd2.str() + "}"; 
  
  const double chisq_dt_b0b_svd1 = get_chisq(h_dt_b0b_svd1, h_dtb0b_pdf_svd1, 3);
  std::stringstream ss_chisq_dt_b0b_svd1;
  ss_chisq_dt_b0b_svd1 << std::setprecision(2) << chisq_dt_b0b_svd1;
  const std::string schisq_dt_b0b_svd1 = "#scale[1.5]{#chi^{2}/ndf=" + ss_chisq_dt_b0b_svd1.str() + "}";
  
  const double chisq_dt_b0b_svd2 = get_chisq(h_dt_b0b_svd2, h_dtb0b_pdf_svd2, 3);
  std::stringstream ss_chisq_dt_b0b_svd2;
  ss_chisq_dt_b0b_svd2 << std::setprecision(2) << chisq_dt_b0b_svd2;
  const std::string schisq_dt_b0b_svd2 = "#scale[1.5]{#chi^{2}/ndf=" + ss_chisq_dt_b0b_svd2.str() + "}"; 
  
  TCanvas *canvas41 = new TCanvas("canvas41","",1024, 768);  
  canvas41->SetLogy(1);
  hist_style(h_dt_b0_svd1, "t_b0_svd1");
  h_dt_b0_svd1->Draw("E1");
  h_dtb0_pdf_svd1->SetLineColor(kBlue);
  TLatex *ltxt_b0_svd1 = new TLatex(sb_dt_ll+0.005, 0.9*h_dt_b0_svd1->GetMaximum(),
			    schisq_dt_b0_svd1.c_str());
  ltxt_b0_svd1->Draw("same");
  h_dtb0_pdf_svd1->Draw("same"); 
  canvas41->Update();
  canvas41->Print("dt_b0_svd1.eps");  
  
  TCanvas *canvas42 = new TCanvas("canvas42","",1024, 768);  
  canvas42->SetLogy(1);
  hist_style(h_dt_b0_svd2, "t_b0_svd2");
  h_dt_b0_svd2->Draw("E1");
  h_dtb0_pdf_svd2->SetLineColor(kBlue);
  TLatex *ltxt_b0_svd2 = new TLatex(sb_dt_ll+0.005, 0.9*h_dt_b0_svd2->GetMaximum(),
			    schisq_dt_b0_svd2.c_str());
  ltxt_b0_svd2->Draw("same");
  h_dtb0_pdf_svd2->Draw("same"); 
  canvas42->Update();
  canvas42->Print("dt_b0_svd2.eps");
   
  TCanvas *canvas51 = new TCanvas("canvas51","",1024, 768);
  canvas51->SetLogy(1);
  hist_style(h_dt_b0b_svd1, "t_b0b_svd1");
  h_dt_b0b_svd1->Draw("E1");
  h_dtb0b_pdf_svd1->Draw("same");
  h_dtb0b_pdf_svd1->SetLineColor(kBlue);
  TLatex *ltxt_b0b_svd1 = new TLatex(sb_dt_ll+0.005, 0.9*h_dt_b0_svd1->GetMaximum(),
			    schisq_dt_b0b_svd1.c_str());
  ltxt_b0b_svd1->Draw("same");
  canvas51->Update();
  canvas51->Print("dt_b0b_svd1.eps");  
  
  TCanvas *canvas52 = new TCanvas("canvas52","",1024, 768);
  canvas52->SetLogy(1);
  hist_style(h_dt_b0b_svd2, "t_b0b_svd2");
  h_dt_b0b_svd2->Draw("E1");
  h_dtb0b_pdf_svd2->Draw("same");
  h_dtb0b_pdf_svd2->SetLineColor(kBlue);
  TLatex *ltxt_b0b_svd2 = new TLatex(sb_dt_ll+0.005, 0.9*h_dt_b0_svd2->GetMaximum(),
			    schisq_dt_b0b_svd2.c_str());
  ltxt_b0b_svd2->Draw("same");
  canvas52->Update();
  canvas52->Print("dt_b0b_svd2.eps");
 
  
 const double chisq_dt_b0_svd1_wkp = get_chisq(h_dt_b0_svd1_wkp, h_dtb0_pdf_svd1_wkp, 3);
  std::stringstream ss_chisq_dt_b0_svd1_wkp;
  ss_chisq_dt_b0_svd1_wkp << std::setprecision(2) << chisq_dt_b0_svd1_wkp;
  const std::string schisq_dt_b0_svd1_wkp = "#scale[1.5]{#chi^{2}/ndf=" + ss_chisq_dt_b0_svd1_wkp.str() + "}";
  
  const double chisq_dt_b0_svd2_wkp = get_chisq(h_dt_b0_svd2_wkp, h_dtb0_pdf_svd2_wkp, 3);
  std::stringstream ss_chisq_dt_b0_svd2_wkp;
  ss_chisq_dt_b0_svd2_wkp << std::setprecision(2) << chisq_dt_b0_svd2_wkp;
  const std::string schisq_dt_b0_svd2_wkp = "#scale[1.5]{#chi^{2}/ndf=" + ss_chisq_dt_b0_svd2_wkp.str() + "}"; 
  
  const double chisq_dt_b0b_svd1_wkp = get_chisq(h_dt_b0b_svd1_wkp, h_dtb0b_pdf_svd1_wkp, 3);
  std::stringstream ss_chisq_dt_b0b_svd1_wkp;
  ss_chisq_dt_b0b_svd1_wkp << std::setprecision(2) << chisq_dt_b0b_svd1_wkp;
  const std::string schisq_dt_b0b_svd1_wkp = "#scale[1.5]{#chi^{2}/ndf=" + ss_chisq_dt_b0b_svd1_wkp.str() + "}";
  
  const double chisq_dt_b0b_svd2_wkp = get_chisq(h_dt_b0b_svd2_wkp, h_dtb0b_pdf_svd2_wkp, 3);
  std::stringstream ss_chisq_dt_b0b_svd2_wkp;
  ss_chisq_dt_b0b_svd2_wkp << std::setprecision(2) << chisq_dt_b0b_svd2_wkp;
  const std::string schisq_dt_b0b_svd2_wkp = "#scale[1.5]{#chi^{2}/ndf=" + ss_chisq_dt_b0b_svd2_wkp.str() + "}"; 
  
  TCanvas *canvas41_wkp = new TCanvas("canvas41_wkp","",1024, 768);  
  canvas41_wkp->SetLogy(1);
  hist_style(h_dt_b0_svd1_wkp, "t_b0_svd1_wkp");
  h_dt_b0_svd1_wkp->Draw("E1");
  h_dtb0_pdf_svd1_wkp->SetLineColor(kBlue);
  TLatex *ltxt_b0_svd1_wkp = new TLatex(sb_dt_ll+0.005, 0.9*h_dt_b0_svd1_wkp->GetMaximum(),
			    schisq_dt_b0_svd1_wkp.c_str());
  ltxt_b0_svd1_wkp->Draw("same");
  h_dtb0_pdf_svd1_wkp->Draw("same"); 
  canvas41_wkp->Update();
  canvas41_wkp->Print("dt_b0_svd1_wkp.eps");  
  
  TCanvas *canvas42_wkp = new TCanvas("canvas42_wkp","",1024, 768);  
  canvas42_wkp->SetLogy(1);
  hist_style(h_dt_b0_svd2_wkp, "t_b0_svd2_wkp");
  h_dt_b0_svd2_wkp->Draw("E1");
  h_dtb0_pdf_svd2_wkp->SetLineColor(kBlue);
  TLatex *ltxt_b0_svd2_wkp = new TLatex(sb_dt_ll+0.005, 0.9*h_dt_b0_svd2_wkp->GetMaximum(),
			    schisq_dt_b0_svd2_wkp.c_str());
  ltxt_b0_svd2_wkp->Draw("same");
  h_dtb0_pdf_svd2_wkp->Draw("same"); 
  canvas42_wkp->Update();
  canvas42_wkp->Print("dt_b0_svd2_wkp.eps");
   
  TCanvas *canvas51_wkp = new TCanvas("canvas51_wkp","",1024, 768);
  canvas51_wkp->SetLogy(1);
  hist_style(h_dt_b0b_svd1_wkp, "t_b0b_svd1_wkp");
  h_dt_b0b_svd1_wkp->Draw("E1");
  h_dtb0b_pdf_svd1_wkp->Draw("same");
  h_dtb0b_pdf_svd1_wkp->SetLineColor(kBlue);
  TLatex *ltxt_b0b_svd1_wkp = new TLatex(sb_dt_ll+0.005, 0.9*h_dt_b0_svd1_wkp->GetMaximum(),
			    schisq_dt_b0b_svd1_wkp.c_str());
  ltxt_b0b_svd1_wkp->Draw("same_wkp");
  canvas51_wkp->Update();
  canvas51_wkp->Print("dt_b0b_svd1_wkp.eps");  
  
  TCanvas *canvas52_wkp = new TCanvas("canvas52_wkp","",1024, 768);
  canvas52_wkp->SetLogy(1);
  hist_style(h_dt_b0b_svd2_wkp, "t_b0b_svd2_wkp");
  h_dt_b0b_svd2_wkp->Draw("E1");
  h_dtb0b_pdf_svd2_wkp->Draw("same");
  h_dtb0b_pdf_svd2_wkp->SetLineColor(kBlue);
  TLatex *ltxt_b0b_svd2_wkp = new TLatex(sb_dt_ll+0.005, 0.9*h_dt_b0_svd2_wkp->GetMaximum(),
			    schisq_dt_b0b_svd2_wkp.c_str());
  ltxt_b0b_svd2_wkp->Draw("same");
  canvas52_wkp->Update();
  canvas52_wkp->Print("dt_b0b_svd2_wkp.eps");  
  
  /**/
  const double chisq_fd_svd1 = get_chisq(h_fd_svd1, h_fd_pdf_svd1, 5);
  std::stringstream ss_chisq_fd_svd1;
  ss_chisq_fd_svd1 << std::setprecision(2) << chisq_fd_svd1;
  const std::string schisq_fd_svd1 = "#scale[1.5]{#chi^{2}/ndf=" + ss_chisq_fd_svd1.str() + "}";
  
  const double chisq_fd_svd2 = get_chisq(h_fd_svd2, h_fd_pdf_svd2, 5);
  std::stringstream ss_chisq_fd_svd2;
  ss_chisq_fd_svd2 << std::setprecision(2) << chisq_fd_svd2;
  const std::string schisq_fd_svd2 = "#scale[1.5]{#chi^{2}/ndf=" + ss_chisq_fd_svd2.str() + "}";  
  
  TCanvas *canvas61 = new TCanvas("canvas61","",1024, 768);
  normres_style(canvas61);
  canvas61->cd(1);
  h_fd_svd1->Draw("E1");
  h_fd_pdf_svd1->SetLineColor(kBlue);
  h_fd_pdf_svd1->Draw("C SAME");
  TLatex *ltxt_fd_svd1 = new TLatex(h_fd_ll+0.005, 0.9*h_fd_svd1->GetMaximum(),
			    schisq_fd_svd1.c_str());
  ltxt_fd_svd1->Draw();
  canvas61->cd(2);
  normres_style(r_fd1, normres_sigma, "F_{B#bar{B}/q#bar{q}}", "");
  r_fd1->Draw("E1");
  for( unsigned int i=0; i<normres_sigma.size(); i++ )
    normres_sigma.at(i).Draw(); 
  canvas61->Update();
  canvas61->Print("fd_svd1.eps");
  
  TCanvas *canvas62 = new TCanvas("canvas62","",1024, 768);
  normres_style(canvas62);
  canvas62->cd(1);
  h_fd_svd2->Draw("E1");
  h_fd_pdf_svd2->SetLineColor(kBlue);
  h_fd_pdf_svd2->Draw("C SAME");
  TLatex *ltxt_fd_svd2 = new TLatex(h_fd_ll+0.005, 0.9*h_fd_svd2->GetMaximum(),
			    schisq_fd_svd2.c_str());
  ltxt_fd_svd2->Draw();
  canvas62->cd(2);
  normres_style(r_fd2, normres_sigma, "F_{B#bar{B}/q#bar{q}}", "");
  r_fd2->Draw("E1");
  for( unsigned int i=0; i<normres_sigma.size(); i++ )
    normres_sigma.at(i).Draw(); 
  canvas62->Update();
  canvas62->Print("fd_svd2.eps"); 
  
  
   TCanvas *canvas66 = new TCanvas("canvas66","",1024, 768);
   canvas66->Divide(4,2);
   TCanvas *canvas67 = new TCanvas("canvas67","",1024, 768);
   canvas67->Divide(4,2); 
   for( unsigned int i=0; i<bins_rbin; i++ )
     plot_histogram(fcn, par, canvas66, canvas67, i, data10, "wks");
 
   canvas66->Update();
   canvas66->Print("fd_rbin_svd1.eps");
   canvas67->Update();
   canvas67->Print("fd_rbin_svd2.eps"); 

  const double chisq_fd_svd1_wkp = get_chisq(h_fd_svd1_wkp, h_fd_pdf_svd1_wkp, 5);
  std::stringstream ss_chisq_fd_svd1_wkp;
  ss_chisq_fd_svd1_wkp << std::setprecision(2) << chisq_fd_svd1_wkp;
  const std::string schisq_fd_svd1_wkp = "#scale[1.5]{#chi^{2}/ndf=" + ss_chisq_fd_svd1_wkp.str() + "}";
  
  const double chisq_fd_svd2_wkp = get_chisq(h_fd_svd2_wkp, h_fd_pdf_svd2_wkp, 5);
  std::stringstream ss_chisq_fd_svd2_wkp;
  ss_chisq_fd_svd2_wkp << std::setprecision(2) << chisq_fd_svd2_wkp;
  const std::string schisq_fd_svd2_wkp = "#scale[1.5]{#chi^{2}/ndf=" + ss_chisq_fd_svd2_wkp.str() + "}";  
  
  TCanvas *canvas61_wkp = new TCanvas("canvas61_wkp","",1024, 768);
  normres_style(canvas61_wkp);
  canvas61_wkp->cd(1);
  h_fd_svd1_wkp->Draw("E1");
  h_fd_pdf_svd1_wkp->SetLineColor(kBlue);
  h_fd_pdf_svd1_wkp->Draw("C SAME");
  TLatex *ltxt_fd_svd1_wkp = new TLatex(h_fd_ll+0.005, 0.9*h_fd_svd1_wkp->GetMaximum(),
			    schisq_fd_svd1_wkp.c_str());
  ltxt_fd_svd1_wkp->Draw();
  canvas61_wkp->cd(2);
  normres_style(r_fd1_wkp, normres_sigma, "F_{B#bar{B}/q#bar{q}}", "");
  r_fd1_wkp->Draw("E1");
  for( unsigned int i=0; i<normres_sigma.size(); i++ )
    normres_sigma.at(i).Draw(); 
  canvas61_wkp->Update();
  canvas61_wkp->Print("fd_svd1_wkp.eps");
  
  TCanvas *canvas62_wkp = new TCanvas("canvas62_wkp","",1024, 768);
  normres_style(canvas62_wkp);
  canvas62_wkp->cd(1);
  h_fd_svd2_wkp->Draw("E1");
  h_fd_pdf_svd2_wkp->SetLineColor(kBlue);
  h_fd_pdf_svd2_wkp->Draw("C SAME");
  TLatex *ltxt_fd_svd2_wkp = new TLatex(h_fd_ll+0.005, 0.9*h_fd_svd2_wkp->GetMaximum(),
			    schisq_fd_svd2_wkp.c_str());
  ltxt_fd_svd2_wkp->Draw();
  canvas62_wkp->cd(2);
  normres_style(r_fd2_wkp, normres_sigma, "F_{B#bar{B}/q#bar{q}}", "");
  r_fd2_wkp->Draw("E1");
  for( unsigned int i=0; i<normres_sigma.size(); i++ )
    normres_sigma.at(i).Draw(); 
  canvas62_wkp->Update();
  canvas62_wkp->Print("fd_svd2_wkp.eps"); 
  
  
   TCanvas *canvas66_wkp = new TCanvas("canvas66_wkp","",1024, 768);
   canvas66_wkp->Divide(4,2);
   TCanvas *canvas67_wkp = new TCanvas("canvas67_wkp","",1024, 768);
   canvas67_wkp->Divide(4,2); 
   for( unsigned int i=0; i<bins_rbin; i++ )
     plot_histogram(fcn, par, canvas66_wkp, canvas67_wkp, i, data11, "wkp");
 
   canvas66_wkp->Update();
   canvas66_wkp->Print("fd_rbin_svd1_wkp.eps");
   canvas67_wkp->Update();
   canvas67_wkp->Print("fd_rbin_svd2_wkp.eps");   
   
   
  const double chisq_mbc_svd1 = get_chisq(h_mbc_svd1, h_mbc_svd1_pdf, 8);
  std::stringstream ss_chisq_mbc_svd1;
  ss_chisq_mbc_svd1 << std::setprecision(2) << chisq_mbc_svd1;
  const std::string schisq_mbc_svd1 = "#scale[1.5]{#chi^{2}/ndf=" + ss_chisq_mbc_svd1.str() + "}";
  
  const double chisq_mbc_svd2 = get_chisq(h_mbc_svd2, h_mbc_svd2_pdf, 8);
  std::stringstream ss_chisq_mbc_svd2;
  ss_chisq_mbc_svd2 << std::setprecision(2) << chisq_mbc_svd2;
  const std::string schisq_mbc_svd2 = "#scale[1.5]{#chi^{2}/ndf=" + ss_chisq_mbc_svd2.str() + "}";  

  
  TCanvas *canvas71 = new TCanvas("canvas71","",1024, 768);
  normres_style(canvas71);
  canvas71->cd(1);
  h_mbc_svd1->Draw("E1");
  h_mbc_svd1_pdf->SetLineColor(kBlue);
  h_mbc_svd1_pdf->Draw("C SAME");
  TLatex *ltxt_mbc_svd1 = new TLatex(h_mbc_ll+0.005, 0.9*h_mbc_svd1->GetMaximum(),
			    schisq_mbc_svd1.c_str());
  ltxt_mbc_svd1->Draw();  
  canvas71->cd(2);
  normres_style(r_mbc1, normres_sigma, "B^{0} M_{bc}", "GeV/c^{2}");
  r_mbc1->Draw("E1");
  for( unsigned int i=0; i<normres_sigma.size(); i++ )
    normres_sigma.at(i).Draw(); 
  canvas71->Update();
  canvas71->Print("mbc_svd1_wks.eps");
  
  TCanvas *canvas72 = new TCanvas("canvas72","",1024, 768);
  normres_style(canvas72);
  canvas72->cd(1);
  h_mbc_svd2->Draw("E1");
  h_mbc_svd2_pdf->SetLineColor(kBlue);
  h_mbc_svd2_pdf->Draw("C SAME");
  TLatex *ltxt_mbc_svd2 = new TLatex(h_mbc_ll+0.005, 0.9*h_mbc_svd2->GetMaximum(),
			    schisq_mbc_svd2.c_str());
  ltxt_mbc_svd2->Draw();  
  canvas72->cd(2);
  normres_style(r_mbc2, normres_sigma, "B^{0} M_{bc}", "GeV/c^{2}");
  r_mbc2->Draw("E1");
  for( unsigned int i=0; i<normres_sigma.size(); i++ )
    normres_sigma.at(i).Draw(); 
  canvas72->Update();
  canvas72->Print("mbc_svd2_wks.eps");

  const double chisq_mbc_svd1_wkp = get_chisq(h_mbc_svd1_wkp, h_mbc_svd1_pdf_wkp, 8);
  std::stringstream ss_chisq_mbc_svd1_wkp;
  ss_chisq_mbc_svd1_wkp << std::setprecision(2) << chisq_mbc_svd1_wkp;
  const std::string schisq_mbc_svd1_wkp = "#scale[1.5]{#chi^{2}/ndf=" + ss_chisq_mbc_svd1_wkp.str() + "}";
  
  const double chisq_mbc_svd2_wkp = get_chisq(h_mbc_svd2_wkp, h_mbc_svd2_pdf_wkp, 8);
  std::stringstream ss_chisq_mbc_svd2_wkp;
  ss_chisq_mbc_svd2_wkp << std::setprecision(2) << chisq_mbc_svd2_wkp;
  const std::string schisq_mbc_svd2_wkp= "#scale[1.5]{#chi^{2}/ndf=" + ss_chisq_mbc_svd2_wkp.str() + "}";  
  
  TCanvas *canvas71_wkp = new TCanvas("canvas71_wkp","",1024, 768);
  normres_style(canvas71_wkp);
  canvas71_wkp->cd(1);
  h_mbc_svd1_wkp->Draw("E1");
  h_mbc_svd1_pdf_wkp->SetLineColor(kBlue);
  h_mbc_svd1_pdf_wkp->Draw("C SAME");
  TLatex *ltxt_mbc_svd1_wkp = new TLatex(h_mbc_ll+0.005, 0.9*h_mbc_svd1_wkp->GetMaximum(),
			    schisq_mbc_svd1_wkp.c_str());
  ltxt_mbc_svd1_wkp->Draw();  
  canvas71_wkp->cd(2);
  normres_style(r_mbc1_wkp, normres_sigma, "B^{+} M_{bc}", "GeV/c^{2}");
  r_mbc1_wkp->Draw("E1");
  for( unsigned int i=0; i<normres_sigma.size(); i++ )
    normres_sigma.at(i).Draw(); 
  canvas71_wkp->Update();
  canvas71_wkp->Print("mbc_svd1_wkp.eps");
  
  TCanvas *canvas72_wkp = new TCanvas("canvas72_wkp","",1024, 768);
  normres_style(canvas72_wkp);
  canvas72_wkp->cd(1);
  h_mbc_svd2_wkp->Draw("E1");
  h_mbc_svd2_pdf_wkp->SetLineColor(kBlue);
  h_mbc_svd2_pdf_wkp->Draw("C SAME");
  TLatex *ltxt_mbc_svd2_wkp = new TLatex(h_mbc_ll+0.005, 0.9*h_mbc_svd2_wkp->GetMaximum(),
			    schisq_mbc_svd2_wkp.c_str());
  ltxt_mbc_svd2_wkp->Draw();  
  canvas72_wkp->cd(2);
  normres_style(r_mbc2_wkp, normres_sigma, "B^{+} M_{bc}", "GeV/c^{2}");
  r_mbc2_wkp->Draw("E1");
  for( unsigned int i=0; i<normres_sigma.size(); i++ )
    normres_sigma.at(i).Draw(); 
  canvas72_wkp->Update();
  canvas72_wkp->Print("mbc_svd2_wkp.eps");
  
  TCanvas *canvas431 = new TCanvas("canvas431","",1024, 1024);
  normres_style(canvas431);
  canvas431->cd(1); 
  hist_style(h_de_mw_proj1, "m_{3#pi}", "GeV/c^{2}");
  h_de_mw_proj1->SetMinimum(0.01);
  h_de_mw_proj1->Draw("E1");
  // h_de_mw_proj1->SetLineColor(kBlue);
  p_de_mw_proj1->Draw("HIST SAME");
  p_de_mw_proj1->SetLineColor(kBlue);
  // p_de_mw_proj1->SetLineStyle(kDashed);

  canvas431->cd(2);
  normres_style(r_de1, normres_sigma, "m_{3#pi}", "GeV/c^{2}");
  r_de1->Draw("E1");
  for( unsigned int i=0; i<normres_sigma.size(); i++ )
    normres_sigma.at(i).Draw();

  canvas431->Update();
  canvas431->Print("de_mw_proj_1_wks.eps");   
   
  TCanvas *canvas432 = new TCanvas("canvas432","",1024, 1024);
  normres_style(canvas432);
  canvas432->cd(1);

  hist_style(h_de_mw_proj2, "M_{3#pi}", "GeV/c^{2}");
  h_de_mw_proj2->SetMinimum(0.01);
 h_de_mw_proj2->Draw("E1");
  //  h_de_mw_proj2->SetLineColor(kBlue);
  p_de_mw_proj2->Draw("HIST SAME");
  p_de_mw_proj2->SetLineColor(kBlue);
  // p_de_mw_proj2->SetLineStyle(kDashed);

  canvas432->cd(2);
  normres_style(r_de2, normres_sigma, "m_{3#pi}", "GeV/c^{2}");
  r_de2->Draw("E1");
  for( unsigned int i=0; i<normres_sigma.size(); i++ )
    normres_sigma.at(i).Draw();


  canvas432->Update();
  canvas432->Print("de_mw_proj_2_wks.eps");   
  
  TCanvas *canvas433 = new TCanvas("canvas433","",1024, 1024);
  normres_style(canvas433);
  canvas433->cd(1);
  h_de_mw_proj3->SetMinimum(0.01);

  hist_style(h_de_mw_proj3, "M_{3#pi}", "GeV/c^{2}");
  h_de_mw_proj3->Draw("E1");
  //  h_de_mw_proj3->SetLineColor(kBlue);
  p_de_mw_proj3->Draw("HIST SAME");
  p_de_mw_proj3->SetLineColor(kBlue);
  //p_de_mw_proj3->SetLineStyle(kDashed);

  canvas433->cd(2);
  normres_style(r_de3, normres_sigma, "m_{3#pi}", "GeV/c^{2}");
  r_de3->Draw("E1");
  for( unsigned int i=0; i<normres_sigma.size(); i++ )
    normres_sigma.at(i).Draw();


  canvas433->Update();
  canvas433->Print("de_mw_proj_3_wks.eps");     
  
  
  TCanvas *canvas434 = new TCanvas("canvas434","",1024, 1024);
  normres_style(canvas434);
  canvas434->cd(1);
  h_de_mw_proj4->SetMinimum(0.01);
hist_style(h_de_mw_proj4, "M_{3#pi}", "GeV/c^{2}");
  h_de_mw_proj4->Draw("E1");
  //  h_de_mw_proj4->SetLineColor(kBlue);
  p_de_mw_proj4->Draw("HIST SAME");
  p_de_mw_proj4->SetLineColor(kBlue);
  // p_de_mw_proj4->SetLineStyle(kDashed);

//    canvas14->cd(2);
//    normres_style(r_hw_sb, normres_sigma, "H_{3#pi}");
//    r_hw_sb->Draw("E1");
//    for( unsigned int i=0; i<normres_sigma.size(); i++ )
//      normres_sigma.at(i).Draw();

  canvas434->cd(2);
  normres_style(r_de4, normres_sigma, "m_{3#pi}", "GeV/c^{2}");
  r_de4->Draw("E1");
  for( unsigned int i=0; i<normres_sigma.size(); i++ )
    normres_sigma.at(i).Draw();


  canvas434->Update();
  canvas434->Print("de_mw_proj_4_wks.eps");     

  TCanvas *canvas435 = new TCanvas("canvas435","",1024, 1024);
//   normres_style(canvas43);
  normres_style(canvas435);
  canvas435->cd(1);
  h_de_mw_proj5->SetMinimum(0.01);
  hist_style(h_de_mw_proj5, "M_{3#pi}", "GeV/c^{2}");
  h_de_mw_proj5->Draw("E1");
  // h_de_mw_proj5->SetLineColor(kBlue);
  p_de_mw_proj5->Draw("HIST SAME");
   p_de_mw_proj5->SetLineColor(kBlue); 
   // p_de_mw_proj5->SetLineStyle(kDashed);

   canvas435->cd(2);
   normres_style(r_de5, normres_sigma, "m_{3#pi}", "GeV/c^{2}");
   r_de5->Draw("E1");
   for( unsigned int i=0; i<normres_sigma.size(); i++ )
     normres_sigma.at(i).Draw();

  canvas435->Update();
  canvas435->Print("de_mw_proj_5_wks.eps");     

  TCanvas *canvas431_wkp = new TCanvas("canvas431_wkp","",1024, 1024);
  normres_style(canvas431_wkp);
  canvas431_wkp->cd(1); 
  hist_style(h_de_mw_proj1_wkp, "m_{3#pi}", "GeV/c^{2}");
  h_de_mw_proj1_wkp->SetMinimum(0.01);
  h_de_mw_proj1_wkp->Draw("E1");
  // h_de_mw_proj1->SetLineColor(kBlue);
  p_de_mw_proj1_wkp->Draw("HIST SAME");
  p_de_mw_proj1_wkp->SetLineColor(kBlue);
  // p_de_mw_proj1->SetLineStyle(kDashed);

  canvas431_wkp->cd(2);
  normres_style(r_de1_wkp, normres_sigma, "m_{3#pi}", "GeV/c^{2}");
  r_de1_wkp->Draw("E1");
  for( unsigned int i=0; i<normres_sigma.size(); i++ )
    normres_sigma.at(i).Draw();

  canvas431_wkp->Update();
  canvas431_wkp->Print("de_mw_proj_1_wkp.eps");   
   
  TCanvas *canvas432_wkp = new TCanvas("canvas432_wkp","",1024, 1024);
  normres_style(canvas432_wkp);
  canvas432_wkp->cd(1);

  hist_style(h_de_mw_proj2_wkp, "M_{3#pi}", "GeV/c^{2}");
  h_de_mw_proj2_wkp->SetMinimum(0.01);
 h_de_mw_proj2_wkp->Draw("E1");
  //  h_de_mw_proj2->SetLineColor(kBlue);
  p_de_mw_proj2_wkp->Draw("HIST SAME");
  p_de_mw_proj2_wkp->SetLineColor(kBlue);
  // p_de_mw_proj2->SetLineStyle(kDashed);

  canvas432_wkp->cd(2);
  normres_style(r_de2_wkp, normres_sigma, "m_{3#pi}", "GeV/c^{2}");
  r_de2_wkp->Draw("E1");
  for( unsigned int i=0; i<normres_sigma.size(); i++ )
    normres_sigma.at(i).Draw();
  canvas432_wkp->Update();
  canvas432_wkp->Print("de_mw_proj_2_wkp.eps");   
  
  TCanvas *canvas433_wkp = new TCanvas("canvas433_wkp","",1024, 1024);
  normres_style(canvas433_wkp);
  canvas433_wkp->cd(1);
  h_de_mw_proj3_wkp->SetMinimum(0.01);

  hist_style(h_de_mw_proj3_wkp, "M_{3#pi}", "GeV/c^{2}");
  h_de_mw_proj3_wkp->Draw("E1");
  //  h_de_mw_proj3->SetLineColor(kBlue);
  p_de_mw_proj3_wkp->Draw("HIST SAME");
  p_de_mw_proj3_wkp->SetLineColor(kBlue);
  //p_de_mw_proj3->SetLineStyle(kDashed);

  canvas433_wkp->cd(2);
  normres_style(r_de3_wkp, normres_sigma, "m_{3#pi}", "GeV/c^{2}");
  r_de3_wkp->Draw("E1");
  for( unsigned int i=0; i<normres_sigma.size(); i++ )
    normres_sigma.at(i).Draw();


  canvas433_wkp->Update();
  canvas433_wkp->Print("de_mw_proj_3_wkp.eps");     
  
  
  TCanvas *canvas434_wkp = new TCanvas("canvas434_wkp","",1024, 1024);
  normres_style(canvas434_wkp);
  canvas434_wkp->cd(1);
  h_de_mw_proj4_wkp->SetMinimum(0.01);
hist_style(h_de_mw_proj4_wkp, "M_{3#pi}", "GeV/c^{2}");
  h_de_mw_proj4_wkp->Draw("E1");
  //  h_de_mw_proj4->SetLineColor(kBlue);
  p_de_mw_proj4_wkp->Draw("HIST SAME");
  p_de_mw_proj4_wkp->SetLineColor(kBlue);
  // p_de_mw_proj4->SetLineStyle(kDashed);

//    canvas14->cd(2);
//    normres_style(r_hw_sb, normres_sigma, "H_{3#pi}");
//    r_hw_sb->Draw("E1");
//    for( unsigned int i=0; i<normres_sigma.size(); i++ )
//      normres_sigma.at(i).Draw();

  canvas434_wkp->cd(2);
  normres_style(r_de4_wkp, normres_sigma, "m_{3#pi}", "GeV/c^{2}");
  r_de4_wkp->Draw("E1");
  for( unsigned int i=0; i<normres_sigma.size(); i++ )
    normres_sigma.at(i).Draw();


  canvas434_wkp->Update();
  canvas434_wkp->Print("de_mw_proj_4_wkp.eps");     

  TCanvas *canvas435_wkp = new TCanvas("canvas435_wkp","",1024, 1024);
//   normres_style(canvas43);
  normres_style(canvas435_wkp);
  canvas435_wkp->cd(1);
  h_de_mw_proj5_wkp->SetMinimum(0.01);
  hist_style(h_de_mw_proj5_wkp, "M_{3#pi}", "GeV/c^{2}");
  h_de_mw_proj5_wkp->Draw("E1");
  // h_de_mw_proj5->SetLineColor(kBlue);
  p_de_mw_proj5_wkp->Draw("HIST SAME");
   p_de_mw_proj5_wkp->SetLineColor(kBlue); 
   // p_de_mw_proj5->SetLineStyle(kDashed);

   canvas435_wkp->cd(2);
   normres_style(r_de5_wkp, normres_sigma, "m_{3#pi}", "GeV/c^{2}");
   r_de5_wkp->Draw("E1");
   for( unsigned int i=0; i<normres_sigma.size(); i++ )
     normres_sigma.at(i).Draw();

  canvas435_wkp->Update();
  canvas435_wkp->Print("de_mw_proj_5_wkp.eps");     
  
  
  return 0;
}

void plot_histogram( Sig3DFcn& fcn, const std::vector<double>& par,
             TCanvas *canvas1, TCanvas *canvas2, const int& rbin, std::vector<Data> data, const std::string& mode )
{
  //Plot histogram
  TH1D *h_fd_svd1 = new TH1D( Form("h_fd_svd1%d",rbin), "", 2*bins_fd, h_fd_ll, h_fd_ul );
  TH1D *h_fd_svd2 = new TH1D( Form("h_fd_svd2%d",rbin), "", 2*bins_fd, h_fd_ll, h_fd_ul );
  TH1D *h_fd_pdf_svd1 =
    new TH1D( Form("h_fd_pdf_svd1_%d",rbin), "", 2*bins_fd, h_fd_ll, h_fd_ul );
  TH1D *h_fd_pdf_svd2 =
    new TH1D( Form("h_fd_pdf_svd2_%d",rbin), "", 2*bins_fd, h_fd_ll, h_fd_ul );
    
if(mode=="wks") {
  for( std::vector<Data>::const_iterator it = fcn.GetData1().begin();
       it != fcn.GetData1().end(); it++ )
    {
      Data data_ = *(it);

      const int r_bin = set_rbin( data_.GetR() );

      if( r_bin == rbin ){
	if( data_.GetExpNo() < 29 )
	  h_fd_svd1->Fill( data_.GetFD() );
	else
	  h_fd_svd2->Fill( data_.GetFD() );
      }
    }
} else {
  for( std::vector<Data>::const_iterator it = fcn.GetData2().begin();
       it != fcn.GetData2().end(); it++ )
    {
      Data data_ = *(it);

      const int r_bin = set_rbin( data_.GetR() );

      if( r_bin == rbin ){
	if( data_.GetExpNo() < 29 )
	  h_fd_svd1->Fill( data_.GetFD() );
	else
	  h_fd_svd2->Fill( data_.GetFD() );
      }
    }
}

if(mode=="wks") {
  for( unsigned int i=0; i < 2*bins_fd; i++ )
    {
      double h_pdf_fd_svd1 = 0.0;
      double h_pdf_fd_svd2 = 0.0;
      
      const double h_fd_svd1 = h_fd_pdf_svd1->GetBinCenter(i+1);
      const double h_fd_svd2 = h_fd_pdf_svd1->GetBinCenter(i+1);
    for( unsigned int l=0; l<data.size(); l++ ){
     int r_bin = set_rbin( data[l].GetR() );
     if(r_bin == rbin){
      if( data[l].GetExpNo() < 29 )
	h_pdf_fd_svd1 += fcn.wksFD(data[l], h_fd_svd1, rbin, par)/fcn.NormwksFD("svd1", rbin, par);
      else
	h_pdf_fd_svd2 += fcn.wksFD(data[l], h_fd_svd2, rbin, par)/fcn.NormwksFD("svd2", rbin, par);    
     }
    }
    h_fd_pdf_svd1->Fill(h_fd_svd1, h_pdf_fd_svd1);
    h_fd_pdf_svd2->Fill(h_fd_svd2, h_pdf_fd_svd2);
    }
    
  h_fd_pdf_svd1->Scale(h_fd_pdf_svd1->GetBinWidth(0));
  h_fd_pdf_svd2->Scale(h_fd_pdf_svd2->GetBinWidth(0));

  std::cout << "rbin svd1 " << rbin << ": chisq = "
<< get_chisq(h_fd_svd1, h_fd_pdf_svd1, 6) << std::endl;  
std::cout << "rbin svd2 " << rbin << ": chisq = "
<< get_chisq(h_fd_svd2, h_fd_pdf_svd2, 6) << std::endl;

  canvas1->cd(rbin+1);

  hist_style(h_fd_svd1, "F_{B#bar{B}/q#bar{q}}");
  h_fd_svd1->Draw("E1");
  h_fd_pdf_svd1->SetLineColor(kBlue);
  h_fd_pdf_svd1->Draw("HIST SAME");

  canvas2->cd(rbin+1);

  hist_style(h_fd_svd2, "F_{B#bar{B}/q#bar{q}}");
  h_fd_svd2->Draw("E1");
  h_fd_pdf_svd2->SetLineColor(kBlue);
  h_fd_pdf_svd2->Draw("HIST SAME");
  
  h_fd_pdf_svd1->Delete();
  h_fd_pdf_svd2->Delete();
}else{
  for( unsigned int i=0; i < 2*bins_fd; i++ )
    {
      double h_pdf_fd_svd1 = 0.0;
      double h_pdf_fd_svd2 = 0.0;
      
      const double h_fd_svd1 = h_fd_pdf_svd1->GetBinCenter(i+1);
      const double h_fd_svd2 = h_fd_pdf_svd1->GetBinCenter(i+1);
    for( unsigned int l=0; l<data.size(); l++ ){
     int r_bin = set_rbin( data[l].GetR() );
     if(r_bin == rbin){
      if( data[l].GetExpNo() < 29 )
	h_pdf_fd_svd1 += fcn.wkpFD(data[l], h_fd_svd1, rbin, par)/fcn.NormwkpFD("svd1", rbin, par);
      else
	h_pdf_fd_svd2 += fcn.wkpFD(data[l], h_fd_svd2, rbin, par)/fcn.NormwkpFD("svd2", rbin, par);    
     }
    }
    h_fd_pdf_svd1->Fill(h_fd_svd1, h_pdf_fd_svd1);
    h_fd_pdf_svd2->Fill(h_fd_svd2, h_pdf_fd_svd2);
    }
    
  h_fd_pdf_svd1->Scale(h_fd_pdf_svd1->GetBinWidth(0));
  h_fd_pdf_svd2->Scale(h_fd_pdf_svd2->GetBinWidth(0));

  std::cout << "rbin svd1 " << rbin << ": chisq = "
<< get_chisq(h_fd_svd1, h_fd_pdf_svd1, 6) << std::endl;  
std::cout << "rbin svd2 " << rbin << ": chisq = "
<< get_chisq(h_fd_svd2, h_fd_pdf_svd2, 6) << std::endl;

  canvas1->cd(rbin+1);

  hist_style(h_fd_svd1, "F_{B#bar{B}/q#bar{q}}");
  h_fd_svd1->Draw("E1");
  h_fd_pdf_svd1->SetLineColor(kBlue);
  h_fd_pdf_svd1->Draw("HIST SAME");

  canvas2->cd(rbin+1);

  hist_style(h_fd_svd2, "F_{B#bar{B}/q#bar{q}}");
  h_fd_svd2->Draw("E1");
  h_fd_pdf_svd2->SetLineColor(kBlue);
  h_fd_pdf_svd2->Draw("HIST SAME");
  
  h_fd_pdf_svd1->Delete();
  h_fd_pdf_svd2->Delete();
}
} 