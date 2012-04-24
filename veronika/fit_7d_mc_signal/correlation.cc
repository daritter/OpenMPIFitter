#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

#include "belle.h"

#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnPrint.h"

#include "TROOT.h"
#include "TCanvas.h"
#include "TH2.h"

#include "tatami/tatami.h"

#include "../Data.h"
#include "../hist.h"
#include "../wtag.h"
#include "../constant.h"

using namespace Belle;

int main( int argc, char* argv[] )
{
  if( argc != 2 )
    {
      std::cout << "ERROR: Invalid number of arguments" << std::endl;
      std::cout << "USAGE: ./correlation type" << std::endl;
      exit(1);
    }

  const std::string type = argv[1];

  if( type != "sig" && type != "mis" && type != "qqbar" &&
      type != "gb0b0b" && type != "gbpbm" &&
      type != "rb0b0b" && type != "rbpbm" &&
      type != "veto" && type != "data" )
    {
      std::cout << "ERROR: Invalid argument supplied" << std::endl;
      std::cout << "USAGE: ./correlation type" << std::endl;
      std::cout << "Valid types are: sig, mis, qqbar, gb0b0b, gbpbm, rb0b0b, rbpbm, a2pi, rho0rho0, b1pi, rho0pipi, pipipipi, data"
		<< std::endl;
      exit(1);
    }

  // Load data
  std::ifstream datafile1, datafile2, datafile3;
  std::ifstream datafile6, datafile7, datafile8;

  if( type == "sig" || type == "mis" )
    {
      datafile1.open("../svd_mc_signal_wkp.txt");
      datafile2.open("../svd_mc_signal_thrust_wkp.txt");
    }
  if( type == "qqbar" )
    {
      datafile1.open("../svd_on_resonance_sideband.txt");
      datafile2.open("../svd_on_resonance_sideband_thrust.txt");
    }
  if( type == "gb0b0b" )
    {
      datafile1.open("../svd_generic_mixed.txt");
      datafile2.open("../svd_generic_mixed_thrust.txt");
    }
  if( type == "gbpbm" )
    {
      datafile1.open("../svd_generic_charged.txt");
      datafile2.open("../svd_generic_charged_thrust.txt");
    }
  if( type == "rb0b0b" )
    {
      datafile1.open("../svd_rareb_mixed.txt");
      datafile2.open("../svd_rareb_mixed_thrust.txt");
    }
  if( type == "rbpbm" )
    {
      datafile1.open("../svd_rareb_charged.txt");
      datafile2.open("../svd_rareb_charged_thrust.txt");
    }
  if( type == "veto" )
    {
      datafile1.open("../svd_generic_vetoed_modes.txt");
      datafile2.open("../svd_generic_vetoed_modes_thrust.txt");
    }
  if( type == "data" )
    {
      datafile1.open("../svd_on_resonance.txt");
      datafile2.open("../svd_on_resonance_thrust.txt");
    }

  std::vector<Data> data1;
  if( type == "sig" )
    {
      FillPlot(datafile1, datafile2, data1, 1, 0, 1);
    }
  else if( type == "mis" )
    {
      FillPlot(datafile1, datafile2, data1, 1, 0, 2);
    }
  else if( type == "qqbar" )
    {
      FillPlot(datafile1, datafile2, data1, 1, 0, 0, ctbto_cut, side_mbc_ll, side_mbc_ul, side_de_ll, side_de_ul, h_momega_ll, h_momega_ul);
    }
  else
    {
      FillPlot(datafile1, datafile2, data1, 1, 0, 0);
    }

  //Global Settings
  set_style();

  double de_ul = h_de_ul, de_ll = h_de_ll;

  //Create histogram
  if( type == "qqbar" ){
    de_ul = side_de_ul, de_ll = side_de_ll;
  }

  // double  h_dt_ll = -7.5, h_dt_ul = 7.5; 
  
  TH2D *h_defd    = new TH2D( "h_defd", "", bins_de, de_ll, de_ul,
				bins_fd, h_fd_ll, h_fd_ul );
  TH2D *h_demw   = new TH2D( "h_demw", "", bins_de, de_ll, de_ul,
			      bins_momega, h_momega_ll, h_momega_ul );
  TH2D *h_dehw   = new TH2D( "h_dehw", "", bins_de, de_ll, de_ul,
				bins_homega, h_homega_ll, h_homega_ul );  
  TH2D *h_dedt   = new TH2D( "h_dedt", "", bins_de, de_ll, de_ul,
				 bins_dt, h_dt_ll, h_dt_ul );  
  TH2D *h_derbin  = new TH2D( "h_derbin", "", bins_de, de_ll, de_ul,
				bins_rbin, rbin_bound );  
  TH2D *h_dembc  = new TH2D( "h_dembc", "", bins_de, de_ll, de_ul,
				 bins_mbc, h_mbc_ll, h_mbc_ul ); 
  TH2D *h_fdmw   = new TH2D( "h_fdmw", "", bins_fd, h_fd_ll, h_fd_ul,
			      bins_momega, h_momega_ll, h_momega_ul );
  TH2D *h_fdhw   = new TH2D( "h_fdhw", "", bins_fd, h_fd_ll, h_fd_ul,
			      bins_homega, h_homega_ll, h_homega_ul );
  TH2D *h_fddt  = new TH2D( "h_fddt", "", bins_fd, h_fd_ll, h_fd_ul,
			      bins_dt, h_dt_ll, h_dt_ul );  
  TH2D *h_fdrbin  = new TH2D( "h_fdrbin", "", bins_fd, h_fd_ll, h_fd_ul,
			      bins_rbin, rbin_bound );
  TH2D *h_fdmbc  = new TH2D( "h_fdmbc", "", bins_fd, h_fd_ll, h_fd_ul,
			      bins_mbc, h_mbc_ll, h_mbc_ul ); 
  TH2D *h_mwhw  = new TH2D( "h_mwhw", "", bins_momega, h_momega_ll, h_momega_ul,
			      bins_homega, h_homega_ll, h_homega_ul );
  TH2D *h_mwdt   = new TH2D( "h_mwdt", "", bins_de, h_momega_ll, h_momega_ul,
			      bins_dt, h_dt_ll, h_dt_ul ); 
  TH2D *h_mwrbin = new TH2D( "h_mwrbin", "", bins_momega, h_momega_ll, h_momega_ul,
			      bins_rbin, rbin_bound );
  TH2D *h_mwmbc  = new TH2D( "h_mwmbc", "", bins_momega, h_momega_ll, h_momega_ul,
			      bins_mbc, h_mbc_ll, h_mbc_ul ); 
  TH2D *h_hwdt   = new TH2D( "h_hwdt", "", bins_homega, h_homega_ll, h_homega_ul,
			      bins_dt, h_dt_ll, h_dt_ul );
  TH2D *h_hwrbin = new TH2D( "h_hwrbin", "", bins_homega, h_homega_ll, h_homega_ul,
			      bins_rbin, rbin_bound );
  TH2D *h_hwmbc  = new TH2D( "h_hwmbc", "", bins_homega, h_homega_ll, h_homega_ul,
			      bins_mbc, h_mbc_ll, h_mbc_ul ); 
  TH2D *h_dtrbin   = new TH2D( "h_dtrbin", "", bins_dt, h_dt_ll, h_dt_ul,
			      bins_rbin, rbin_bound   );  
  TH2D *h_dtmbc   = new TH2D( "h_dtmbc", "", bins_dt, h_dt_ll, h_dt_ul,
			      bins_mbc, h_mbc_ll, h_mbc_ul  );  
  TH2D *h_mbcrbin   = new TH2D( "h_mbcrbin", "", bins_mbc, h_mbc_ll, h_mbc_ul,
			      bins_rbin, rbin_bound   );  
			      
  for( std::vector<Data>::iterator it = data1.begin();
       it != data1.end(); it++ )
    {
      Data data_ = *(it);

      h_defd->Fill( data_.GetDE(), data_.GetFD() );
      h_demw->Fill( data_.GetDE(), data_.GetMomega() );
      h_dehw->Fill( data_.GetDE(), data_.GetHomega() );
      h_derbin->Fill( data_.GetDE(), data_.GetR() );
      h_dembc->Fill( data_.GetDE(), data_.GetMbc() );
      h_fdmw->Fill( data_.GetFD(), data_.GetMomega() );
      h_fdhw->Fill( data_.GetFD(), data_.GetHomega() );
      h_fdrbin->Fill( data_.GetFD(), data_.GetR() );
      h_fdmbc->Fill( data_.GetFD(), data_.GetMbc() );
      h_mwhw->Fill( data_.GetMomega(), data_.GetHomega() );
      h_mwrbin->Fill( data_.GetMomega(), data_.GetR() );
      h_mwmbc->Fill( data_.GetMomega(), data_.GetMbc() );
      h_hwrbin->Fill( data_.GetHomega(), data_.GetR() );
      h_hwmbc->Fill( data_.GetHomega(), data_.GetMbc() );
      h_dedt->Fill( data_.GetDE(), data_.GetDeltat() );
      h_fddt->Fill( data_.GetFD(), data_.GetDeltat() );
      h_mwdt->Fill( data_.GetMomega(), data_.GetDeltat() );
      h_hwdt->Fill( data_.GetHomega(), data_.GetDeltat() );
      h_dtmbc->Fill( data_.GetDeltat(), data_.GetMbc() );
      h_dtrbin->Fill( data_.GetDeltat(), data_.GetR() );
      h_mbcrbin->Fill( data_.GetMbc(), data_.GetR() );
    }

  std::cout << setiosflags(std::ios::left) << std::setprecision(2)
	    << std::setw(10) << ""
	    << std::setw(10) << "DE"
	    << std::setw(10) << "FD"
	    << std::setw(10) << "Mw"
	    << std::setw(10) << "Hw"
	    << std::setw(10) << "Dt"
	    << std::setw(10) << "R-bin" 
	    << std::setw(10) << "Mbc"	    
	    << std::endl;
  std::cout << setiosflags(std::ios::left) << std::setprecision(2)
	    << std::setw(10) << "DE"
	    << std::setw(10) << 1.0
	    << std::setw(10) << h_defd->GetCorrelationFactor()
	    << std::setw(10) << h_demw->GetCorrelationFactor()
	    << std::setw(10) << h_dehw->GetCorrelationFactor()
	    << std::setw(10) << h_dedt->GetCorrelationFactor()
	    << std::setw(10) << h_derbin->GetCorrelationFactor() 
	    << std::setw(10) << h_dembc->GetCorrelationFactor() 
	    << std::endl;
  std::cout << setiosflags(std::ios::left) << std::setprecision(2)
	    << std::setw(10) << "FD"
	    << std::setw(10) << ""
	    << std::setw(10) << 1.0
	    << std::setw(10) << h_fdmw->GetCorrelationFactor()
	    << std::setw(10) << h_fdhw->GetCorrelationFactor()
	    << std::setw(10) << h_fddt->GetCorrelationFactor()
	    << std::setw(10) << h_fdrbin->GetCorrelationFactor() 
	    << std::setw(10) << h_fdmbc->GetCorrelationFactor() 
	    << std::endl;
  std::cout << setiosflags(std::ios::left) << std::setprecision(2)
	    << std::setw(10) << "Mw"
	    << std::setw(10) << ""
	    << std::setw(10) << ""
	    << std::setw(10) << 1.0
	    << std::setw(10) << h_mwhw->GetCorrelationFactor()
	    << std::setw(10) << h_mwdt->GetCorrelationFactor()
            << std::setw(10) << h_mwrbin->GetCorrelationFactor() 
	    << std::setw(10) << h_mwmbc->GetCorrelationFactor() 
            << std::endl;
  std::cout << setiosflags(std::ios::left) << std::setprecision(2)
	    << std::setw(10) << "Hw"
	    << std::setw(10) << ""
	    << std::setw(10) << ""
	    << std::setw(10) << ""
	    << std::setw(10) << 1.0
	    << std::setw(10) << h_hwdt->GetCorrelationFactor()
	    << std::setw(10) << h_hwrbin->GetCorrelationFactor() 	   
	    << std::setw(10) << h_hwmbc->GetCorrelationFactor() 
	    << std::endl;
  std::cout << setiosflags(std::ios::left) << std::setprecision(2)
	    << std::setw(10) << "Dt"
	    << std::setw(10) << ""
	    << std::setw(10) << ""
	    << std::setw(10) << ""
	    << std::setw(10) << ""
	    << std::setw(10) << 1.0
	    << std::setw(10) << h_dtrbin->GetCorrelationFactor() 
	    << std::setw(10) << h_dtmbc->GetCorrelationFactor() 
	    << std::endl;
  std::cout << setiosflags(std::ios::left) << std::setprecision(2)
	    << std::setw(10) << "R-bin"
	    << std::setw(10) << ""
	    << std::setw(10) << ""
	    << std::setw(10) << ""
	    << std::setw(10) << ""	    
 	    << std::setw(10) << ""
	    << std::setw(10) << 1.0 
	    << std::setw(10) << h_mbcrbin->GetCorrelationFactor() 
	    << std::endl;
  std::cout << setiosflags(std::ios::left) << std::setprecision(2)
	    << std::setw(10) << "Mbc"
	    << std::setw(10) << ""
	    << std::setw(10) << ""
	    << std::setw(10) << ""
	    << std::setw(10) << ""	    
 	    << std::setw(10) << ""
	    << std::setw(10) << "" 
	    << std::setw(10) << 1.0
	    << std::endl;

  TCanvas *canvas1 = new TCanvas("canvas1","",1050, 768);
  canvas1->Divide(2,3);

  TH1D* h_de = h_defd->ProjectionX();
  TH1D* h_fd = h_defd->ProjectionY();
  TH1D* h_mw = h_demw->ProjectionY();
  TH1D* h_hw = h_dehw->ProjectionY();
  TH1D* h_dt = h_dedt->ProjectionY();

  canvas1->cd(1);
  hist_style(h_de, "#DeltaE", "GeV");
  h_de->Draw("E1");
  canvas1->cd(2);
  hist_style(h_fd, "F");
  h_fd->Draw("E1");
  canvas1->cd(3);
  hist_style(h_mw, "m_{3#pi}", "GeV/c^{2}");
  h_mw->Draw("E1");
  canvas1->cd(4);
  hist_style(h_hw, "H_{3#pi}");
  h_hw->Draw("E1");
  canvas1->cd(5);
  hist_style(h_hw, "#Delta t");
  h_dt->Draw("E1");
  
  canvas1->Update();
  canvas1->Print("wks.eps");

  TCanvas *canvas2 = new TCanvas("canvas2","",1200, 768);
  hist_style(h_demw, "#DeltaE [GeV]");
  h_demw->GetYaxis()->SetTitle("m_{3#pi} [GeV/c^{2}]");
  h_demw->Draw("COLZ");
  canvas2->Update();
  canvas2->Print("de_mw_corr.eps");

  TCanvas *canvas3 = new TCanvas("canvas3","",1200, 768);
  hist_style(h_defd, "#DeltaE [GeV]");
  h_defd->GetYaxis()->SetTitle("F_{B#bar{B}/q#bar{q}}");
  h_defd->Draw("COLZ");
  canvas3->Update();
  canvas3->Print("de_fd_corr.eps");

  TCanvas *canvas4 = new TCanvas("canvas4","",1200, 768);
  hist_style(h_dedt, "#DeltaE [GeV]");
  h_defd->GetYaxis()->SetTitle("#Delta t [ps]");
  h_defd->Draw("COLZ");
  canvas4->Update();
  canvas4->Print("de_dt_corr.eps");

  TCanvas *canvas5 = new TCanvas("canvas5","",1200, 768);
  hist_style(h_hwdt, "H_{3#pi}");
  h_hwdt->GetYaxis()->SetTitle("#Delta t [ps]");
  h_hwdt->Draw("COLZ");
  canvas5->Update();
  canvas5->Print("hw_dt_corr.eps");


  return 0;
}
