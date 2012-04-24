#ifndef MN_Sig3DFcn_H_
#define MN_Sig3DFcn_H_

#include "belle.h"
#include "tatami/tatami.h"

#include "Minuit2/FCNBase.h"

using namespace ROOT::Minuit2;

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

namespace ROOT {

namespace Minuit2 {

  class Sig3DFcn : public FCNBase {

  public:
      /** Define how different pieces of the pdf should be reduced to one value */
      typedef std::plus<double> operator_type;
      /** Define the increase of the objective function value which is eqvivalent to 1 sigma */
      static const double error_def = 1.0;

  Sig3DFcn( std::vector<Data>& data1_, std::vector<Data>& data2_ ) :

    data1_(data1_), data2_(data2_),

    //1sigma errors
    error_def_(1.0) {}

  ~Sig3DFcn() {}

  virtual double Up() const { return error_def_; }
  virtual double operator()( const std::vector<double>& parameters ) const;

  /** finalize the pdf after all values from the different processes are collected */
  double finalize(const std::vector<double>& parameters, double value) const;
  /** Load the chunk of data to be used by this process of the pdf
    * Each process will get its own copy of the pdf. Once this is done, the
    * data should be loaded by this method. The responsibilty to divide the
    * data into distinct portions lies with the pdf and should be done here.
    * The first parameter is the process id and the second parameter is the
    * total number of processes.
    *
    * @param process the id of this process, starting from 0
    * @param size the number of processes in total
    */
  void load(int process, int size);

   //Accessors
  void FillYield( const std::string& svd_no, const std::string& mode,
		  std::vector<double>& vNwks,
		  const std::vector<double>& par ) const;

  double wksPDF( const Data& data,
		 const double& de, const double& fd,
		 const double& wm, const double& mbc, const double& wh,
		 const double& dt, const int& q,
		 const int& r_bin,
		 const std::vector<double>& par ) const;
  std::vector<double> NormwksPDF( const std::string& svd_no,
				  const std::vector<double>& par ) const;

  double wksDE( const Data& data,
		const double& de,
		const std::vector<double>& par,
		unsigned int par0=0 ) const;
  double NormwksDE( const std::string& svd_no,
		    const std::vector<double>& par,
		    unsigned int par0=0 ) const;
  double wksFD( const Data& data,
		const double& fd, const int& r_bin,
		const std::vector<double>& par,
		unsigned int par0=12 ) const;
  double NormwksFD( const std::string& svd_no, const int& r_bin,
		    const std::vector<double>& par,
		    unsigned int par0=12 ) const;
  double wksMw( const Data& data,const std::string& mode,
		const double& wm, const double& de,
		const std::vector<double>& par,
		unsigned int par0=68 ) const;
  double NormwksMw( const std::string& svd_no, const std::string& mode, const double& de,
		    const std::vector<double>& par,
		    unsigned int par0=68 ) const;
  double wksHw( const Data& data,
		const double& wh,
		const std::vector<double>& par,
		unsigned int par0=78 ) const;
  double NormwksHw( const std::string& svd_no,
		    const std::vector<double>& par,
		    unsigned int par0=78 ) const;
  double wksMBC( const Data& data,
		const double& mbc,
		const std::vector<double>& par,
		unsigned int par0=82 ) const;
  double NormwksMBC( const std::string& svd_no,
		    const std::vector<double>& par,
		    unsigned int par0=82 ) const;
  double wksDtQ( const Data& data,
		 const double& dt, const int& q,
		 const std::vector<double>& par,
		 unsigned int par0=90 ) const;
  void AddOutDtQ( const Data& data, const double& dt,
		  double& pdf ) const;

  double wkpPDF( const Data& data,
		 const double& de, const double& fd,
		 const double& wm, const double& mbc, const double& wh,
		 const double& dt, const int& q,
		 const int& r_bin,
		 const std::vector<double>& par ) const;
  std::vector<double> NormwkpPDF( const std::string& svd_no,
				  const std::vector<double>& par ) const;

  double wkpDE( const Data& data,
		const double& de,
		const std::vector<double>& par,
		unsigned int par0=0 ) const;//was:208 when not sharing parameters
  double NormwkpDE( const std::string& svd_no,
		    const std::vector<double>& par,
		    unsigned int par0=0 ) const;//was:208 when not sharing parameters
  double wkpFD( const Data& data,
		const double& fd, const int& r_bin,
		const std::vector<double>& par,
		unsigned int par0=12 ) const;//was:220 when not sharing parameters
  double NormwkpFD( const std::string& svd_no, const int& r_bin,
		    const std::vector<double>& par,
		    unsigned int par0=12 ) const;//was:220 when not sharing parameters
  double wkpMw( const Data& data,
		const double& wm, const double& de,
		const std::vector<double>& par,
		unsigned int par0=276 ) const;
  double NormwkpMw( const std::string& svd_no, const double& de,
		    const std::vector<double>& par,
		    unsigned int par0=276 ) const;
  double wkpHw( const Data& data,
		const double& wh,
		const std::vector<double>& par,
		unsigned int par0=286 ) const;
  double NormwkpHw( const std::string& svd_no,
		    const std::vector<double>& par,
		    unsigned int par0=286 ) const;
  double wkpMBC( const Data& data,
		const double& mbc,
		const std::vector<double>& par,
		unsigned int par0=82 ) const;//was:288 when not sharing parameters
  double NormwkpMBC( const std::string& svd_no,
		    const std::vector<double>& par,
		    unsigned int par0=82 ) const;//was:288 when not sharing parameters
  double wkpDtQ( const Data& data,
		 const double& dt, const int& q,
		 const std::vector<double>& par,
		 unsigned int par0=302 ) const;
  void PlotMwPDF( const double& de, const double& mw,
			  double& mw_sig_pdf, const std::vector<double>& par);
  void PlotMwPDF_wkp( const double& de, const double& mw,
			  double& mw_sig_pdf, const std::vector<double>& par);

  unsigned int GetNpar() const {return npar_wks_svd1_;}
  const std::vector<Data>& GetData1() const {return data1_;}
  const std::vector<Data>& GetData2() const {return data2_;}
/*
  //Modifiers
  void SetXll( double de_ll ) {de_ll_ = de_ll;}
  void SetXul( double de_ul ) {de_ul_ = de_ul;}
  void SetYll( double om_ll ) {om_ll_ = om_ll;}
  void SetYul( double om_ul ) {om_ul_ = om_ul;}
  void SetZll( double oh_ll ) {oh_ll_ = oh_ll;}
  void SetZul( double oh_ul ) {oh_ul_ = oh_ul;}
  void SetTll( double dt_ll ) {dt_ll_ = dt_ll;}
  void SetTul( double dt_ul ) {dt_ul_ = dt_ul;}
  void SetQll( double q_ll ) {q_ll_ = q_ll;}
  void SetQul( double q_ul ) {q_ul_ = q_ul;}
  void SetLRll( double fd_ll ) {fd_ll_ = fd_ll;}
  void SetLRul( double fd_ul ) {fd_ul_ = fd_ul;}
*/
  void SetErrorDef( double def ) {error_def_ = def;}

  private:
  std::vector<Data>& data1_;
  std::vector<Data>& data2_;
/*
  double de_ll_;
  double de_ul_;
  double om_ll_;
  double om_ul_;
  double oh_ll_;
  double oh_ul_;
  double dt_ll_;
  double dt_ul_;
  double q_ll_;
  double q_ul_;
  double fd_ll_;
  double fd_ul_;
*/
  static const unsigned int npar_wks_svd1_ = 106;
  static const unsigned int npar_wkp_svd1_ = 106;

  double error_def_;
  };

} // namespace Minuit2

} // namespace ROOT

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

#endif //MN_Sig1DFcn_H_
