#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>

#include "Minuit2/FunctionMinimum.h"
#include "belle.h"
#include "tatami/tatami.h"
#include "../func.h"
#include "../cheb.h"
#include "../tools.h"
#include "../wtag.h"
#include "../hist.h"
#include "Sig3DFcn.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif
using namespace ROOT::Minuit2;


double Sig3DFcn::operator()( const std::vector<double>& par ) const
{
//     std::cout << "par.size(): " << par.size() << std::endl;
 assert( par.size() == (npar_wks_svd1_+npar_wkp_svd1_)*2 );

  //Normalise
  const std::vector<double> int_wks_pdf_svd1 = NormwksPDF("svd1", par);
  const std::vector<double> int_wks_pdf_svd2 = NormwksPDF("svd2", par);
  const std::vector<double> int_wkp_pdf_svd1 = NormwkpPDF("svd1", par);
  const std::vector<double> int_wkp_pdf_svd2 = NormwkpPDF("svd2", par);

  //Yields
  std::vector<double> Nwks_svd1, Nwks_svd2, Nwkp_svd1, Nwkp_svd2;
  FillYield("svd1", "wks", Nwks_svd1, par);
  FillYield("svd2", "wks", Nwks_svd2, par);
  FillYield("svd1", "wkp", Nwkp_svd1, par);
  FillYield("svd2", "wkp", Nwkp_svd2, par);

  double log_pdf_wks = 0.0, log_pdf_wkp = 0.0, log_pdf = 0.0;
  unsigned int i;
  int r_bin1, r_bin2;
  double wks_pdf, pdf, wkp_pdf;
#pragma omp parallel private(i,data1,r_bin1,wks_pdf,pdf) reduction(+:log_pdf_wks)
  {
#pragma omp for
    for( i=0; i<data1_.size(); i++ )
      {
	Data &data1 = data1_[i];
	//R-bin
	r_bin1 = set_rbin( data1.GetR() );

	//PI+ PI- component
	wks_pdf = wksPDF( data1, data1.GetDE(), data1.GetFD(), data1.GetMomega(), data1.GetMbc(), data1.GetHomega(), data1.GetDeltat(), data1.GetQ(), r_bin1, par );

	//PDF
	if( data1.GetExpNo() < 29 ){
	  pdf = (Nwks_svd1[r_bin1]*wks_pdf/int_wks_pdf_svd1[r_bin1]);
	}else{
	  pdf = (Nwks_svd2[r_bin1]*wks_pdf/int_wks_pdf_svd2[r_bin1]);
	}
	if( 0.0 < pdf )
	  log_pdf_wks += log( pdf );
      }
  }
#pragma omp parallel private(i,data2,r_bin2,wkp_pdf,pdf) reduction(+:log_pdf_wkp)
  {
#pragma omp for
    for( i=0; i<data2_.size(); i++ )
      {
	Data &data2 = data2_[i];
	//R-bin
	r_bin2 = set_rbin( data2.GetR() );

	//PI+ PI- component
	wkp_pdf = wkpPDF( data2, data2.GetDE(), data2.GetFD(), data2.GetMomega(), data2.GetMbc(), data2.GetHomega(), data2.GetDeltat(), data2.GetQ(), r_bin2, par );

	//PDF
	if( data2.GetExpNo() < 29 ){
	  pdf = Nwkp_svd1[r_bin2] *wkp_pdf/int_wkp_pdf_svd1[r_bin2];
	    //std::cout << "Nwkp_svd1[r_bin2] " << Nwkp_svd1[r_bin2]  << std::endl;
	}else{
	  pdf = Nwkp_svd2[r_bin2] *wkp_pdf/int_wkp_pdf_svd2[r_bin2];
	    //std::cout << "Nwkp_svd2[r_bin2] " << Nwkp_svd2[r_bin2]  << std::endl;
	}
	if( 0.0 < pdf )
	  log_pdf_wkp += log( pdf );
      }

  }
  return log_pdf_wks + log_pdf_wkp;
}

double Sig3DFcn::finalize( const std::vector<double>& par, double value) const
{
  //Yields
  std::vector<double> Nwks_svd1, Nwks_svd2, Nwkp_svd1, Nwkp_svd2;
  FillYield("svd1", "wks", Nwks_svd1, par);
  FillYield("svd2", "wks", Nwks_svd2, par);
  FillYield("svd1", "wkp", Nwkp_svd1, par);
  FillYield("svd2", "wkp", Nwkp_svd2, par);

  double log_pdf = value - Nwks_svd1[7] - Nwks_svd2[7] - Nwkp_svd1[7] - Nwkp_svd2[7];

  std::cout << std::setprecision(16) << "-2logL = " << -2.0*log_pdf << std::endl;

  return ( -2.0*log_pdf );
}

template<class T> void splitData(std::vector<T> &data, int i, int N){
    const size_t size = std::ceil(1.0 * data.size() / N);
    const size_t start = i * size;
    const size_t end = std::min(data.size(), start + size);
    std::vector<T> tmp;
    tmp.reserve(end-start);
    tmp.insert(tmp.end(),data.begin()+start,data.begin()+end);
    std::swap(tmp,data);
}

void Sig3DFcn::load(int process, int size){
    if(size<=1) return;
    splitData(data1_,process,size);
    splitData(data2_,process,size);
}

void Sig3DFcn::FillYield( const std::string& svd_no, const std::string& mode,
			  std::vector<double>& vNwks,
			  const std::vector<double>& par ) const
{
  unsigned int Nwks_par0 = 0;
  if( mode == "wkp"){
    Nwks_par0 += 2*npar_wks_svd1_;
      if( svd_no == "svd2")
	Nwks_par0 += npar_wkp_svd1_;
  }else{
      if( svd_no == "svd2")
	Nwks_par0 += npar_wks_svd1_;
  }

  int add;

  if( mode == "wkp") add = 92;
  else add = 92;

  const double Nwks  = par[Nwks_par0+add+13];
  const double Nwks0 = Nwks*par[Nwks_par0+add+1]*par[Nwks_par0+add+7];
  const double Nwks1 = Nwks*par[Nwks_par0+add+2]*par[Nwks_par0+add+8];
  const double Nwks2 = Nwks*par[Nwks_par0+add+3]*par[Nwks_par0+add+9];
  const double Nwks3 = Nwks*par[Nwks_par0+add+4]*par[Nwks_par0+add+10];
  const double Nwks4 = Nwks*par[Nwks_par0+add+5]*par[Nwks_par0+add+11];
  const double Nwks5 = Nwks*par[Nwks_par0+add+6]*par[Nwks_par0+add+12];
  const double Nwks6 =
    Nwks*(1.0-(par[Nwks_par0+add+1]*par[Nwks_par0+add+7])-
	   (par[Nwks_par0+add+2]*par[Nwks_par0+add+8])-
	   (par[Nwks_par0+add+3]*par[Nwks_par0+add+9])-
	   (par[Nwks_par0+add+4]*par[Nwks_par0+add+10])-
	   (par[Nwks_par0+add+5]*par[Nwks_par0+add+11])-
	   (par[Nwks_par0+add+6]*par[Nwks_par0+add+12]));

 /*
  std::cout << "MODE " << mode;
  std::cout << " " << svd_no << std::endl;
  std::cout << "Nwks0 " << Nwks0 << std::endl;
  std::cout << "Nwks1 " << Nwks1 << std::endl;
  std::cout << "Nwks2 " << Nwks2 << std::endl;
  std::cout << "Nwks3 " << Nwks3 << std::endl;
  std::cout << "Nwks4 " << Nwks4 << std::endl;
  std::cout << "Nwks5 " << Nwks5 << std::endl;
  std::cout << "Nwks6 " << Nwks6 << std::endl;
  std::cout << "Nwks " << Nwks << std::endl;
*/
  vNwks.push_back(Nwks0);
  vNwks.push_back(Nwks1);
  vNwks.push_back(Nwks2);
  vNwks.push_back(Nwks3);
  vNwks.push_back(Nwks4);
  vNwks.push_back(Nwks5);
  vNwks.push_back(Nwks6);
  vNwks.push_back(Nwks);

  return;
}

double Sig3DFcn::wksPDF( const Data& data,
			const double& de, const double& fd,
			const double& wm, const double& mbc, const double& wh,
			const double& dt, const int& q,
			const int& r_bin,
			const std::vector<double>& par ) const
{
    std::string svd_no;
  if( data.GetExpNo() > 29 )
    svd_no = "svd2";
  else
    svd_no = "svd1";

    const double wks_pdf =
    wksDE(data, de, par)*
    wksFD(data, fd, r_bin, par)*
   ( wksMw(data, "wks", wm, de, par)/NormwksMw(svd_no, "wks", de, par))*
    wksMBC(data, mbc, par)*
    wksDtQ(data, dt, q, par )*
    wksHw(data, wh, par)
    ;

  return wks_pdf;
}

double Sig3DFcn::wkpPDF( const Data& data,
			const double& de, const double& fd,
			const double& wm, const double& mbc, const double& wh,
			const double& dt, const int& q,
			const int& r_bin,
			const std::vector<double>& par ) const
{
    std::string svd_no;
  if( data.GetExpNo() > 29 )
    svd_no = "svd2";
  else
    svd_no = "svd1";

    const double wks_pdf =
    wkpDE(data, de, par)*
    wkpFD(data, fd, r_bin, par)*
   ( wksMw(data, "wkp", wm, de, par)/NormwksMw(svd_no, "wkp", de, par))*
    wkpMBC(data, mbc, par)*
    wkpDtQ(data, dt, q, par )*
   wksHw(data, wh, par)*
   1.0
    ;

  return wks_pdf;
}

std::vector<double>
Sig3DFcn::NormwksPDF( const std::string& svd_no,
		       const std::vector<double>& par ) const
{
  std::vector<double> int_wks_pdf;
  for( unsigned int i=0; i<bins_rbin; i++ )
    int_wks_pdf.push_back( NormwksDE(svd_no, par)*
			   NormwksFD(svd_no, i, par)*
			  // do not uncomment me NormwksMw(svd_no, de, par)*
			   NormwksMBC(svd_no, par)*
			   NormwksHw(svd_no, par)*
			   1.0);

  return int_wks_pdf;
}

std::vector<double>
Sig3DFcn::NormwkpPDF( const std::string& svd_no,
		       const std::vector<double>& par ) const
{
  std::vector<double> int_wks_pdf;
  for( unsigned int i=0; i<bins_rbin; i++ )
    int_wks_pdf.push_back( NormwkpDE(svd_no, par)*
			   NormwkpFD(svd_no, i, par)*
			   // do not uncomment me NormwksMw(svd_no, de, par)*
			  NormwkpMBC(svd_no, par)*
			   NormwksHw(svd_no, par)*
			   1.0);

  return int_wks_pdf;
}

double Sig3DFcn::wksDE( const Data& data,
		       const double& de,
		       const std::vector<double>& par,
		       unsigned int par0 ) const
{
  if( data.GetExpNo() > 29 )
    par0 += npar_wks_svd1_;

  const double mean1  = par[par0+0];
  const double sigma1 = par[par0+1];
  const double mean2  = par[par0+2] + mean1;
  const double sigma2 = par[par0+3] * sigma1;
  const double mean3  = par[par0+4] + mean2;
  const double sigma3 = par[par0+5] * sigma2;
  const double c = par[par0+6];
  const double frac1  = par[par0+7];
  const double frac2  = par[par0+8];
  const double frac3  = par[par0+9];

  double wks_pdf = 0;

    wks_pdf =
    frac1*gaussian( de, mean1, sigma1 ) +
    frac2*gaussian( de, mean2, sigma2 )
    +(1. - frac1 - frac2 - frac3)*gaussian( de, mean3, sigma3 )
    +frac3*chebyshev1(de, c);

  return wks_pdf;
}

double Sig3DFcn::wkpDE( const Data& data,
		       const double& de,
		       const std::vector<double>& par,
		       unsigned int par0 ) const
{
  if( data.GetExpNo() > 29 )
    par0 += 2*npar_wks_svd1_+npar_wkp_svd1_;
  else{
    par0 += 2*npar_wks_svd1_;
  }
  const double mean1  = par[par0+0];
  const double sigma1 = par[par0+1];

  par0 -= 2*npar_wks_svd1_;

  const double mean2  = par[par0+2] + mean1;
  const double sigma2 = par[par0+3] * sigma1;
  const double mean3  = par[par0+4] + mean2;
  const double sigma3 = par[par0+5] * sigma2;
  const double c = par[par0+6];
  const double frac1  = par[par0+7];
  const double frac2  = par[par0+8];
  const double frac3  = par[par0+9];

  double wks_pdf = 0;

  wks_pdf =
    frac1*gaussian( de, mean1, sigma1 ) +
    frac2*gaussian( de, mean2, sigma2 )
    +(1. - frac1 - frac2 - frac3)*gaussian( de, mean3, sigma3 )
    +frac3*chebyshev1(de, c);

  return wks_pdf;
}

double Sig3DFcn::NormwksDE( const std::string& svd_no,
			   const std::vector<double>& par,
			   unsigned int par0 ) const
{
  if( svd_no == "svd2" )
    par0 += npar_wks_svd1_;

  const double mean1  = par[par0+0];
  const double sigma1 = par[par0+1];
  const double mean2  = par[par0+2] + mean1;
  const double sigma2 = par[par0+3] * sigma1;
  const double mean3  = par[par0+4] + mean2;
  const double sigma3 = par[par0+5] * sigma2;
  const double c = par[par0+6];
  const double frac1  = par[par0+7];
  const double frac2  = par[par0+8];
  const double frac3  = par[par0+9];

  double int_wks_pdf = 0;

    int_wks_pdf =
    frac1*norm_gaussian( h_de_ll, h_de_ul, mean1, sigma1 )
    + frac2*norm_gaussian( h_de_ll, h_de_ul, mean2, sigma2 )
    + (1. - frac1 - frac2 - frac3)*norm_gaussian( h_de_ll, h_de_ul, mean3, sigma3 )
    +frac3*norm_chebyshev1(h_de_ll, h_de_ul, c);

  return int_wks_pdf;
}

double Sig3DFcn::NormwkpDE( const std::string& svd_no,
			   const std::vector<double>& par,
			   unsigned int par0 ) const
{
  if(  svd_no == "svd2"  )
    par0 += 2*npar_wks_svd1_+npar_wkp_svd1_;
  else{
    par0 += 2*npar_wks_svd1_;
  }
  const double mean1  = par[par0+0];
  const double sigma1 = par[par0+1];

  par0 -= 2*npar_wks_svd1_;

  const double mean2  = par[par0+2] + mean1;
  const double sigma2 = par[par0+3] * sigma1;
  const double mean3  = par[par0+4] + mean2;
  const double sigma3 = par[par0+5] * sigma2;
  const double c = par[par0+6];
  const double frac1  = par[par0+7];
  const double frac2  = par[par0+8];
  const double frac3  = par[par0+9];

  double int_wks_pdf = 0;

    int_wks_pdf =
    frac1*norm_gaussian( h_de_ll, h_de_ul, mean1, sigma1 )
    + frac2*norm_gaussian( h_de_ll, h_de_ul, mean2, sigma2 )
    + (1. - frac1 - frac2 - frac3)*norm_gaussian( h_de_ll, h_de_ul, mean3, sigma3 )
    +frac3*norm_chebyshev1(h_de_ll, h_de_ul, c);



  return int_wks_pdf;
}

double Sig3DFcn::wksFD( const Data& data,
		       const double& fd, const int& r_bin,
		       const std::vector<double>& par,
		       unsigned int par0 ) const
{
  if( data.GetExpNo() > 29 )
    par0 += npar_wks_svd1_;

    const double mean1   = par[par0+(8*r_bin)+0];
    const double sigma1  = par[par0+(8*r_bin)+1];
    const double mean2   = mean1 + par[par0+(8*r_bin)+2];
    const double sigma2  = sigma1 * par[par0+(8*r_bin)+3];
    const double mean3   = mean2 + par[par0+(8*r_bin)+4];
    const double sigma3  = sigma2 * par[par0+(8*r_bin)+5];
    const double frac1   = par[par0+(8*r_bin)+6];
    const double frac2   = par[par0+(8*r_bin)+7];

    const double wks_pdf =
      frac1*gaussian( fd, mean1, sigma1 ) +
      frac2*gaussian( fd, mean2, sigma2 )
      +(1. - frac1 - frac2)*gaussian( fd, mean3, sigma3 );

    return wks_pdf;
}

double Sig3DFcn::NormwksFD( const std::string& svd_no,
			   const int& r_bin,
			   const std::vector<double>& par,
			   unsigned int par0 ) const
{
  if( svd_no == "svd2" )
    par0 += npar_wks_svd1_;
   const double mean1   = par[par0+(8*r_bin)+0];
  const double sigma1  = par[par0+(8*r_bin)+1];
  const double mean2   = mean1 + par[par0+(8*r_bin)+2];
  const double sigma2  = sigma1 * par[par0+(8*r_bin)+3];
  const double mean3   = mean2 + par[par0+(8*r_bin)+4];
  const double sigma3  = sigma2 * par[par0+(8*r_bin)+5];
  const double frac1   = par[par0+(8*r_bin)+6];
  const double frac2   = par[par0+(8*r_bin)+7];

  const double int_wks_pdf =
    frac1*norm_gaussian( h_fd_ll, h_fd_ul, mean1, sigma1 )
    + frac2*norm_gaussian( h_fd_ll, h_fd_ul, mean2, sigma2 )
    + (1. - frac1 - frac2)*norm_gaussian( h_fd_ll, h_fd_ul, mean3, sigma3 );
  return int_wks_pdf;
}

double Sig3DFcn::wkpFD( const Data& data,
		       const double& fd, const int& r_bin,
		       const std::vector<double>& par,
		       unsigned int par0 ) const
{
  if( data.GetExpNo() > 29 )
    par0 += 2*npar_wks_svd1_+npar_wkp_svd1_;
  else{
    par0 += 2*npar_wks_svd1_;
  }
  const double mean1  = par[par0+(8*r_bin)+0];
  const double sigma1 = par[par0+(8*r_bin)+1];

  par0 -= 2*npar_wks_svd1_;

    const double mean2   = mean1 + par[par0+(8*r_bin)+2];
    const double sigma2  = sigma1 * par[par0+(8*r_bin)+3];
    const double mean3   = mean2 + par[par0+(8*r_bin)+4];
    const double sigma3  = sigma2 * par[par0+(8*r_bin)+5];
    const double frac1   = par[par0+(8*r_bin)+6];
    const double frac2   = par[par0+(8*r_bin)+7];

    const double wks_pdf =
      frac1*gaussian( fd, mean1, sigma1 ) +
      frac2*gaussian( fd, mean2, sigma2 )
      +(1. - frac1 - frac2)*gaussian( fd, mean3, sigma3 );

    return wks_pdf;

}

double Sig3DFcn::NormwkpFD( const std::string& svd_no,
			   const int& r_bin,
			   const std::vector<double>& par,
			   unsigned int par0 ) const
{
  if( svd_no == "svd2" )
    par0 += 2*npar_wks_svd1_+npar_wkp_svd1_;
  else{
    par0 += 2*npar_wks_svd1_;
  }
  const double mean1  = par[par0+(8*r_bin)+0];
  const double sigma1 = par[par0+(8*r_bin)+1];

  par0 -= 2*npar_wks_svd1_;

  const double mean2   = mean1 + par[par0+(8*r_bin)+2];
  const double sigma2  = sigma1 * par[par0+(8*r_bin)+3];
  const double mean3   = mean2 + par[par0+(8*r_bin)+4];
  const double sigma3  = sigma2 * par[par0+(8*r_bin)+5];
  const double frac1   = par[par0+(8*r_bin)+6];
  const double frac2   = par[par0+(8*r_bin)+7];

  const double int_wks_pdf =
    frac1*norm_gaussian( h_fd_ll, h_fd_ul, mean1, sigma1 )
    + frac2*norm_gaussian( h_fd_ll, h_fd_ul, mean2, sigma2 )
    + (1. - frac1 - frac2)*norm_gaussian( h_fd_ll, h_fd_ul, mean3, sigma3 );
  return int_wks_pdf;

}


double Sig3DFcn::wksMw( const Data& data,const std::string& mode,
		       const double& wm, const double& de,
		       const std::vector<double>& par,
		       unsigned int par0 ) const{


  if( data.GetExpNo() > 29)
    par0 += npar_wkp_svd1_;

  double mean1 = par[par0+0];
  double sigma1 = par[par0+1];

  if( mode == "wkp"){
    par0 += 2*npar_wks_svd1_;
  }
//  if( data.GetExpNo() > 29 ){
   // par0 += npar_wks_svd1_;
    mean1 +=  par[par0+8]*de;
    sigma1 += par[par0+9]*de*de;
 /* }else{
    mean1 = par[par0+0] + par[par0+8]*de;
    sigma1 = par[par0+1] + par[par0+9]*de*de;
  }*/
  if( mode == "wkp")
    par0 -= 2*npar_wks_svd1_;

  const double mean2  = par[par0+2] + mean1;
  const double sigma2 = par[par0+3] * sigma1;
  const double mean3  = par[par0+4] + mean2;
  const double sigma3 = par[par0+5] * sigma2;
  const double frac1  = par[par0+6];
  const double frac2  = par[par0+7];

   const double wks_pdf =
     frac1*gaussian( wm, mean1, sigma1 ) +
     frac2*gaussian( wm, mean2, sigma2 )
     +(1. - frac1 - frac2)*gaussian( wm, mean3, sigma3 );
  return wks_pdf;
}

double Sig3DFcn::NormwksMw( const std::string& svd_no, const std::string& mode, const double& de,
			   const std::vector<double>& par,
			   unsigned int par0 ) const
{
// std::cout<<"Mode: " << mode << std::endl;
// std::cout<<"SVD: " << svd_no << std::endl;
// //std::cout<<"par0 1: " << par0 << std::endl;
  if( svd_no == "svd2")
    par0 += npar_wks_svd1_;
 //std::cout<<"par0 2: " << par0 << std::endl;

  double mean1 = par[par0+0];
  double sigma1 = par[par0+1];

  if( mode == "wkp"){
    par0 += 2*npar_wks_svd1_;
  }
//std::cout<<"par0 3: " << par0 << std::endl;

  //if( svd_no == "svd2" ){
   // par0 += npar_wks_svd1_;
    mean1 += par[par0+8]*de;
    sigma1 += par[par0+9]*de*de;
 /* }else{
    mean1 = par[par0+0] + par[par0+8]*de;
    sigma1 = par[par0+1] + par[par0+9]*de*de;
  }*/

  if( mode == "wkp")
    par0 -= 2*npar_wks_svd1_;

  const double mean2  = par[par0+2] + mean1;
  const double sigma2 = par[par0+3] * sigma1;
  const double mean3  = par[par0+4] + mean2;
  const double sigma3 = par[par0+5] * sigma2;
  const double frac1  = par[par0+6];
  const double frac2  = par[par0+7];

  const double int_wks_pdf =
    frac1*norm_gaussian( h_momega_ll, h_momega_ul, mean1, sigma1 )
    + frac2*norm_gaussian( h_momega_ll, h_momega_ul, mean2, sigma2 )
    + (1. - frac1 - frac2)*norm_gaussian( h_momega_ll, h_momega_ul, mean3, sigma3 );

  return int_wks_pdf;
}

double Sig3DFcn::wksMBC( const Data& data,
		       const double& mbc,
		       const std::vector<double>& par,
		       unsigned int par0 ) const
{
  if( data.GetExpNo() > 29 )
    par0 += npar_wks_svd1_;

  const double mean1 = par[par0+0];
  const double sigma1 = par[par0+1];
  const double mean2  = par[par0+2] + mean1;
  const double sigma2 = par[par0+3] * sigma1;
  const double mean3  = par[par0+4] + mean2;
  const double sigma3 = par[par0+5] * sigma2;
  const double frac1  = par[par0+6];
  const double frac2  = par[par0+7];

   const double wks_pdf =
     frac1*gaussian( mbc, mean1, sigma1 ) +
     frac2*gaussian( mbc, mean2, sigma2 )
     +(1. - frac1 - frac2)*gaussian( mbc, mean3, sigma3 );

  return wks_pdf;
}

double Sig3DFcn::NormwksMBC( const std::string& svd_no,
			   const std::vector<double>& par,
			   unsigned int par0 ) const
{
  if( svd_no == "svd2" )
    par0 += npar_wks_svd1_;

  const double mean1  = par[par0+0];
  const double sigma1 = par[par0+1];
  const double mean2  = par[par0+2] + mean1;
  const double sigma2 = par[par0+3] * sigma1;
  const double mean3  = par[par0+4] + mean2;
  const double sigma3 = par[par0+5] * sigma2;
  const double frac1  = par[par0+6];
  const double frac2  = par[par0+7];

  const double int_wks_pdf =
    frac1*norm_gaussian( h_mbc_ll, h_mbc_ul, mean1, sigma1 )
    + frac2*norm_gaussian( h_mbc_ll, h_mbc_ul, mean2, sigma2 )
    + (1. - frac1 - frac2)*norm_gaussian( h_mbc_ll, h_mbc_ul, mean3, sigma3 );

  return int_wks_pdf;
}

double Sig3DFcn::wkpMBC( const Data& data,
		       const double& mbc,
		       const std::vector<double>& par,
		       unsigned int par0 ) const
{
  if( data.GetExpNo() > 29 )
    par0 += 2*npar_wks_svd1_+npar_wkp_svd1_;
  else{
    par0 += 2*npar_wks_svd1_;
  }
  const double mean1  = par[par0+0];
  const double sigma1 = par[par0+1];

  par0 -= 2*npar_wks_svd1_;
  const double mean2  = par[par0+2] + mean1;
  const double sigma2 = par[par0+3] * sigma1;
  const double mean3  = par[par0+4] + mean2;
  const double sigma3 = par[par0+5] * sigma2;
  const double frac1  = par[par0+6];
  const double frac2  = par[par0+7];

   const double wks_pdf =
     frac1*gaussian( mbc, mean1, sigma1 ) +
     frac2*gaussian( mbc, mean2, sigma2 )
     +(1. - frac1 - frac2)*gaussian( mbc, mean3, sigma3 );
  return wks_pdf;
}

double Sig3DFcn::NormwkpMBC( const std::string& svd_no,
			   const std::vector<double>& par,
			   unsigned int par0 ) const
{
  if(  svd_no == "svd2"  )
    par0 += 2*npar_wks_svd1_+npar_wkp_svd1_;
  else{
    par0 += 2*npar_wks_svd1_;
  }
  const double mean1  = par[par0+0];
  const double sigma1 = par[par0+1];

  par0 -= 2*npar_wks_svd1_;

  const double mean2  = par[par0+2] + mean1;
  const double sigma2 = par[par0+3] * sigma1;
  const double mean3  = par[par0+4] + mean2;
  const double sigma3 = par[par0+5] * sigma2;
  const double frac1  = par[par0+6];
  const double frac2  = par[par0+7];

  const double int_wks_pdf =
    frac1*norm_gaussian( h_mbc_ll, h_mbc_ul, mean1, sigma1 )
    + frac2*norm_gaussian( h_mbc_ll, h_mbc_ul, mean2, sigma2 )
    + (1. - frac1 - frac2)*norm_gaussian( h_mbc_ll, h_mbc_ul, mean3, sigma3 );

  return int_wks_pdf;
}

double Sig3DFcn::wksHw( const Data& data,
		       const double& wh,
		       const std::vector<double>& par,
		       unsigned int par0 ) const
{
  if( data.GetExpNo() > 29 )
    par0 += npar_wks_svd1_;

  //Prepare coefficients
  std::vector<double> c(4);
  for( unsigned int i=0; i<c.size(); i++ )
    c[i] = par[par0+i];

  const double wks_pdf = chebyshev4( wh, c );

  return wks_pdf;
}

double Sig3DFcn::NormwkpHw( const std::string& svd_no,
			   const std::vector<double>& par,
			   unsigned int par0 ) const
{
  if( svd_no == "svd2" )
    par0 += npar_wkp_svd1_;

  //Prepare coefficients
  std::vector<double> c(4);
  for( unsigned int i=0; i<c.size(); i++ )
    c[i] = par[par0+i];

  const double int_wks_pdf = norm_chebyshev4( h_homega_ll, h_homega_ul, c );

  return int_wks_pdf;
}

double Sig3DFcn::wkpHw( const Data& data,
		       const double& wh,
		       const std::vector<double>& par,
		       unsigned int par0 ) const
{
  if( data.GetExpNo() > 29 )
    par0 += npar_wkp_svd1_;

  //Prepare coefficients
  std::vector<double> c(4);
  for( unsigned int i=0; i<c.size(); i++ )
    c[i] = par[par0+i];

  const double wks_pdf = chebyshev4( wh, c );

  return wks_pdf;
}

double Sig3DFcn::NormwksHw( const std::string& svd_no,
			   const std::vector<double>& par,
			   unsigned int par0 ) const
{
  if( svd_no == "svd2" )
    par0 += npar_wks_svd1_;

  //Prepare coefficients
  std::vector<double> c(4);
  for( unsigned int i=0; i<c.size(); i++ )
    c[i] = par[par0+i];

  const double int_wks_pdf = norm_chebyshev4( h_homega_ll, h_homega_ul, c );

  return int_wks_pdf;
}

double Sig3DFcn::wksDtQ( const Data& data,
			const double& dt, const int& q,
			const std::vector<double>& par,
			unsigned int par0 ) const
{
  if( data.GetExpNo() > 29 )
    par0 += npar_wks_svd1_;

  double Acp_wks, Scp_wks;
  if( data.GetExpNo() < 29 )
    {
      Acp_wks = par[par0+npar_wks_svd1_+1];
      Scp_wks = par[par0+npar_wks_svd1_+2];
    }
  else
    {
      Acp_wks = par[par0+1];
      Scp_wks = par[par0+2];
    }

  const dtres_param_t* const dtres_param =
    get_dtres_param( data.GetExpNo(), data.GetMC() );

  double life_pdf(0), int_life_pdf(0), cos_pdf(0), sin_pdf(0);
  if(!data.getDTComponents(par[par0+0],life_pdf, int_life_pdf, cos_pdf, sin_pdf)){


  //Lifetime component
  life_pdf =      EfRkRdetRnp_fullrec( dt, data.GetBtype(),
				par[par0+0], data.GetAk(), data.GetCk(),
				data.GetRecVNtrk(), data.GetRecVZerr(),
				data.GetRecVChisq(), data.GetRecVNdf(),
				data.GetTagVNtrk(), data.GetTagVZerr(),
				data.GetTagVChisq(), data.GetTagVNdf(),
				data.GetTagVIsl(), dtres_param );

  int_life_pdf =	norm_EfRkRdetRnp_fullrec( h_dt_ll, h_dt_ul, data.GetBtype(),
				par[par0+0], data.GetAk(), data.GetCk(),
				data.GetRecVNtrk(), data.GetRecVZerr(),
				data.GetRecVChisq(), data.GetRecVNdf(),
				data.GetTagVNtrk(), data.GetTagVZerr(),
				data.GetTagVChisq(), data.GetTagVNdf(),
				data.GetTagVIsl(), dtres_param );

  //Acp component
  cos_pdf = 0.5 / par[par0+0] *
      MfRkRdetRnp_fullrec( dt, data.GetBtype(), par[par0+0],
			   data.GetDeltam(), data.GetAk(), data.GetCk(),
			   data.GetRecVNtrk(), data.GetRecVZerr(),
			   data.GetRecVChisq(), data.GetRecVNdf(),
			   data.GetTagVNtrk(), data.GetTagVZerr(),
			   data.GetTagVChisq(), data.GetTagVNdf(),
			   data.GetTagVIsl(), dtres_param );

  //Scp component
  sin_pdf = 0.5 / par[par0+0] *
      AfRkRdetRnp_fullrec( dt, data.GetBtype(), par[par0+0],
			   data.GetDeltam(), data.GetAk(), data.GetCk(),
			   data.GetRecVNtrk(), data.GetRecVZerr(),
			   data.GetRecVChisq(), data.GetRecVNdf(),
			   data.GetTagVNtrk(), data.GetTagVZerr(),
			   data.GetTagVChisq(), data.GetTagVNdf(),
			   data.GetTagVIsl(), dtres_param );
    data.setDTComponents(par[par0+0], life_pdf, int_life_pdf, cos_pdf, sin_pdf);
  }

  //Wrong tag fraction
  const int r_bin = set_rbin( data.GetR() );
  const double wt_w = set_wtag( data.GetExpNo(), data.GetMC(), r_bin );
  const double wt_dw = set_dwtag( data.GetExpNo(), data.GetMC(), r_bin );

  //wks PDF
  double wks_pdf =
    (life_pdf * ( 1.0 - (q*wt_dw) )) +
    (( q * (1.0 - (2.0*wt_w)) ) *
     ((Acp_wks*cos_pdf)+(Scp_wks*sin_pdf)));

  //Normalise over q
  wks_pdf /= 2.0*int_life_pdf;

  AddOutDtQ(data, dt, wks_pdf);

  return wks_pdf;
}

double Sig3DFcn::wkpDtQ( const Data& data,
			const double& dt, const int& q,
			const std::vector<double>& par,
			unsigned int par0 ) const
{
  if( data.GetExpNo() > 29 )
    par0 += npar_wks_svd1_;

  double Acp_wks, Scp_wks;
  if( data.GetExpNo() < 29 )
    {
      Acp_wks = par[par0+npar_wks_svd1_+1];
      Scp_wks = par[par0+npar_wks_svd1_+2];
    }
  else
    {
      Acp_wks = par[par0+1];
      Scp_wks = par[par0+2];
    }

  const dtres_param_t* const dtres_param =
    get_dtres_param( data.GetExpNo(), data.GetMC() );

  double life_pdf(0), int_life_pdf(0), cos_pdf(0), sin_pdf(0);
  if(!data.getDTComponents(par[par0+0],life_pdf, int_life_pdf, cos_pdf, sin_pdf)){
  //Lifetime component
        life_pdf =      EfRkRdetRnp_fullrec( dt, data.GetBtype(),
				par[par0+0], data.GetAk(), data.GetCk(),
				data.GetRecVNtrk(), data.GetRecVZerr(),
				data.GetRecVChisq(), data.GetRecVNdf(),
				data.GetTagVNtrk(), data.GetTagVZerr(),
				data.GetTagVChisq(), data.GetTagVNdf(),
				data.GetTagVIsl(), dtres_param );

        int_life_pdf =	norm_EfRkRdetRnp_fullrec( h_dt_ll, h_dt_ul, data.GetBtype(),
				par[par0+0], data.GetAk(), data.GetCk(),
				data.GetRecVNtrk(), data.GetRecVZerr(),
				data.GetRecVChisq(), data.GetRecVNdf(),
				data.GetTagVNtrk(), data.GetTagVZerr(),
				data.GetTagVChisq(), data.GetTagVNdf(),
				data.GetTagVIsl(), dtres_param );

  //Acp component
        cos_pdf = 0.5 / par[par0+0] *
      MfRkRdetRnp_fullrec( dt, data.GetBtype(), par[par0+0],
			   data.GetDeltam(), data.GetAk(), data.GetCk(),
			   data.GetRecVNtrk(), data.GetRecVZerr(),
			   data.GetRecVChisq(), data.GetRecVNdf(),
			   data.GetTagVNtrk(), data.GetTagVZerr(),
			   data.GetTagVChisq(), data.GetTagVNdf(),
			   data.GetTagVIsl(), dtres_param );

  //Scp component
        sin_pdf = 0.5 / par[par0+0] *
      AfRkRdetRnp_fullrec( dt, data.GetBtype(), par[par0+0],
			   data.GetDeltam(), data.GetAk(), data.GetCk(),
			   data.GetRecVNtrk(), data.GetRecVZerr(),
			   data.GetRecVChisq(), data.GetRecVNdf(),
			   data.GetTagVNtrk(), data.GetTagVZerr(),
			   data.GetTagVChisq(), data.GetTagVNdf(),
			   data.GetTagVIsl(), dtres_param );

    data.setDTComponents(par[par0+0], life_pdf, int_life_pdf, cos_pdf, sin_pdf);
  }

  //Wrong tag fraction
  const int r_bin = set_rbin( data.GetR() );
  const double wt_w = set_wtag( data.GetExpNo(), data.GetMC(), r_bin );
  const double wt_dw = set_dwtag( data.GetExpNo(), data.GetMC(), r_bin );

  //wks PDF
  double wks_pdf =
    (life_pdf * ( 1.0 - (q*wt_dw) )) +
    (( q * (1.0 - (2.0*wt_w)) ) *
     ((Acp_wks*cos_pdf)+(Scp_wks*sin_pdf)));

  //Normalise over q
  wks_pdf /= 2.0*int_life_pdf;

  AddOutDtQ(data, dt, wks_pdf);

  return wks_pdf;
}

void Sig3DFcn::AddOutDtQ( const Data& data, const double& dt,
			 double& pdf ) const
{
  const dtres_param_t* const dtres_param =
    get_dtres_param( data.GetExpNo(), data.GetMC() );

  const double out_pdf =
    gaussian( dt, 0.0, dtres_param->sig_ol ) /
    (2.0*norm_gaussian( h_dt_ll, h_dt_ul, 0.0, dtres_param->sig_ol ));

  double f_ol;
  if( 1 < data.GetRecVNtrk() && 1 < data.GetTagVNtrk() )
    f_ol = dtres_param->fol_mul;
  else
    f_ol = dtres_param->fol_sgl;

  pdf = ( (1.0 - f_ol)*pdf ) + ( f_ol*out_pdf );

  return;
}

void Sig3DFcn::PlotMwPDF( const double& de, const double& mw,
			  double& mw_sig_pdf, const std::vector<double>& par)
{

  Data data_svd1, data_svd2;
  data_svd1.SetExpNo(7);
  data_svd2.SetExpNo(31);

	    mw_sig_pdf += par[105]*wksDE( data_svd1, de, par)
	    *wksMw( data_svd1, "wks", mw, de, par)/(NormwksDE( "svd1", par)
	    *NormwksMw( "svd1", "wks", de, par)) +
	    par[105+npar_wks_svd1_]*wksDE( data_svd2, de, par)
	    *wksMw( data_svd2, "wks", mw, de, par)
	    /(NormwksDE( "svd2", par)
	    *NormwksMw( "svd2", "wks", de, par));

}

void Sig3DFcn::PlotMwPDF_wkp( const double& de, const double& mw,
			  double& mw_sig_pdf, const std::vector<double>& par)
{

  Data data_svd1, data_svd2;
  data_svd1.SetExpNo(7);
  data_svd2.SetExpNo(31);

	    mw_sig_pdf += par[317]*wkpDE( data_svd1, de, par)
	    *wksMw( data_svd1, "wkp", mw, de, par)/(NormwkpDE( "svd1", par)
	    *NormwksMw( "svd1", "wkp", de, par)) +
	    par[317+npar_wks_svd1_]*wkpDE( data_svd2, de, par)
	    *wksMw( data_svd2, "wkp", mw, de, par)
	    /(NormwkpDE( "svd2", par)
	    *NormwksMw( "svd2", "wkp", de, par));

}
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif
