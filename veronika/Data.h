#ifndef Data_H_
#define Data_H_

#include <vector>

#include "constant.h"

#if defined(BELLE_NAMESPACE)
namespace Belle {
#endif

class Data
{
 public:
  Data():cache_lifetime(-1) {};
  ~Data() {}

  //Accessors
  unsigned int GetExpNo() const          {return exp_no_;}
  unsigned int GetRunNo() const          {return run_no_;}
  unsigned int GetEvtNo() const          {return evt_no_;}

  double       GetMbc() const            {return mbc_;}
  double       GetBenergy() const        {return benergy_;}
  double       GetDE() const             {return de_;}
  double       GetLR() const             {return lr_;}
  double       GetFD()  const            {return fd_;}
  double       GetMomega() const         {return momega_;}
  double       GetHomega() const         {return homega_;}

  bool         GetGoodMC() const         {return good_mc_;}
  double       GetMomegaGen() const      {return momega_gen_;}

  double       GetL0c() const            {return l0c_;}
  double       GetL2c() const            {return l2c_;}
  double       GetL0n() const            {return l0n_;}
  double       GetL2n() const            {return l2n_;}
  double       GetCTBTO() const          {return ctbto_;}
  double       GetCTBz() const           {return ctbz_;}
  double       GetCBz() const            {return cbz_;}
  double       GetSumPt() const          {return sumpt_;}

  double       Getpi1Ptot()  const       {return pi1_ptot_;}
  double       Getpi2Ptot() const        {return pi2_ptot_;}
  double       Getpi3Ptot()  const       {return pi3_ptot_;}
  double       Getpi4Ptot()  const       {return pi4_ptot_;}

   double       Getpi1Ctheta() const      {return pi1_ctheta_;}
   double       Getpi2Ctheta() const      {return pi2_ctheta_;}
   double       Getpi3Ctheta() const      {return pi3_ctheta_;}
   double       Getpi4Ctheta() const      {return pi4_ctheta_;}

  double       GetSigPDF() const         {return sig_pdf_;}
  double       GetMisPDF() const         {return mis_pdf_;}
  double       GetQQbarPDF()  const      {return qqbar_pdf_;}
  double       GetGB0B0bPDF() const      {return gb0b0b_pdf_;}
  double       GetGBpBmPDF() const       {return gbpbm_pdf_;}
  double       GetRB0B0bPDF() const      {return rb0b0b_pdf_;}
  double       GetRBpBmPDF() const       {return rbpbm_pdf_;}
  double       Geta2piPDF()  const       {return a2pi_pdf_;}
  double       Getrho0rho0PDF() const    {return rho0rho0_pdf_;}
  double       Getb1piPDF() const        {return b1pi_pdf_;}
  double       Getrho0pipiPDF() const    {return rho0pipi_pdf_;}
  double       GetpipipipiPDF() const    {return pipipipi_pdf_;}

  unsigned int GetMC() const             {return mc_;}
  unsigned int GetBtype() const       {return btype_;}
  double       GetTau()  const           {return tau_;}
  double       GetDeltam()  const        {return delta_m_;}

  double       GetDeltat() const         {return delta_t_;}

  double       GetRecVChisq() const      {return recv_chisq_;}
  unsigned int GetRecVNdf()  const       {return recv_ndf_;}
  double       GetRecVZerr()  const      {return recv_zerr_;}
  unsigned int GetRecVNtrk()  const      {return recv_ntrk_;}

  int          GetQ() const              {return q_;}
  double       GetR()  const             {return r_;}
  int          GetC()  const             {return c_;}

  double       GetTagVChisq()  const     {return tagv_chisq_;}
  unsigned int GetTagVNdf() const        {return tagv_ndf_;}
  double       GetTagVZerr() const       {return tagv_zerr_;}
  unsigned int GetTagVNtrk() const       {return tagv_ntrk_;}
  unsigned int GetTagVIsl() const        {return tagv_isl_;}

  double       GetAk() const             {return a_k_;}
  double       GetCk() const             {return c_k_;}
  double       GetCosTheta() const       {return costheta_;}

  double       GetEfRkRdetRnp() const    {return life_pdf_;}
  double       GetNormEfRkRdetRnp() const {return int_life_pdf_;}
  double       GetMfRkRdetRnp() const    {return cos_pdf_;}
  double       GetAfRkRdetRnp() const    {return sin_pdf_;}

  double       GetSigProb() const        {return sig_prob_;}
  double       GetMisProb() const        {return mis_prob_;}
  double       GetQQbarProb()  const     {return qqbar_prob_;}
  double       GetGB0B0bProb() const     {return gb0b0b_prob_;}
  double       GetGBpBmProb() const      {return gbpbm_prob_;}
  double       GetRB0B0bProb() const     {return rb0b0b_prob_;}
  double       GetRBpBmProb() const      {return rbpbm_prob_;}
  double       Geta2piProb() const       {return a2pi_prob_;}
  double       Getrho0rho0Prob()  const  {return rho0rho0_prob_;}
  double       Getb1piProb() const       {return b1pi_prob_;}
  double       Getrho0pipiProb() const   {return rho0pipi_prob_;}
  double       GetpipipipiProb() const   {return pipipipi_prob_;}

  //Modifiers
  void SetExpNo( const unsigned int& exp_no )         {exp_no_ = exp_no;}
  void SetRunNo( const unsigned int& run_no )         {run_no_ = run_no;}
  void SetEvtNo( const unsigned int& evt_no )         {evt_no_ = evt_no;}

  void SetMbc( const double& mbc )                    {mbc_ = mbc;}
  void SetBenergy( const double& benergy )            {benergy_ = benergy;}
  void SetDE( const double& de )                      {de_ = de;}
  void SetLR( const double& lr )                      {lr_ = lr;}
  void SetFD( const double& fd )                      {fd_ = fd;}
  void SetMomega( const double& momega )                    {momega_ = momega;}
  void SetHomega( const double& homega )                    {homega_ = homega;}
  void SetMRHO0( const double& mrho0 )                {mrho0_ = mrho0;}

  void SetGoodMC( const bool& good_mc )               {good_mc_ = good_mc;}
  void SetMomegaGen( const double& momega_gen )             {momega_gen_ = momega_gen;}

  void SetL0c( const double& l0c )                    {l0c_ = l0c;}
  void SetL2c( const double& l2c )                    {l2c_ = l2c;}
  void SetL0n( const double& l0n )                    {l0n_ = l0n;}
  void SetL2n( const double& l2n )                    {l2n_ = l2n;}
  void SetCTBTO( const double& ctbto )                {ctbto_ = ctbto;}
  void SetCTBz( const double& ctbz )                  {ctbz_ = ctbz;}
  void SetCBz( const double& cbz )                    {cbz_ = cbz;}
  void SetSumPt( const double& sumpt )                {sumpt_ = sumpt;}

  void Setpi1Ptot( const double& pi1_ptot )           {pi1_ptot_ = pi1_ptot;}
  void Setpi2Ptot( const double& pi2_ptot )           {pi2_ptot_ = pi2_ptot;}
  void Setpi3Ptot( const double& pi3_ptot )           {pi3_ptot_ = pi3_ptot;}
  void Setpi4Ptot( const double& pi4_ptot )           {pi4_ptot_ = pi4_ptot;}

  void Setpi1Ctheta( const double& pi1_ctheta )
    {pi1_ctheta_ = pi1_ctheta;}
  void Setpi2Ctheta( const double& pi2_ctheta )
    {pi2_ctheta_ = pi2_ctheta;}
  void Setpi3Ctheta( const double& pi3_ctheta )
    {pi3_ctheta_ = pi3_ctheta;}
  void Setpi4Ctheta( const double& pi4_ctheta )
    {pi4_ctheta_ = pi4_ctheta;}

  void SetSigPDF( const double& sig_pdf )             {sig_pdf_ = sig_pdf;}
  void SetMisPDF( const double& mis_pdf )             {mis_pdf_ = mis_pdf;}
  void SetQQbarPDF( const double& qqbar_pdf )         {qqbar_pdf_ = qqbar_pdf;}
  void SetGB0B0bPDF( const double& gb0b0b_pdf )
    {gb0b0b_pdf_ = gb0b0b_pdf;}
  void SetGBpBmPDF( const double& gbpbm_pdf )         {gbpbm_pdf_ = gbpbm_pdf;}
  void SetRB0B0bPDF( const double& rb0b0b_pdf )
    {rb0b0b_pdf_ = rb0b0b_pdf;}
  void SetRBpBmPDF( const double& rbpbm_pdf )         {rbpbm_pdf_ = rbpbm_pdf;}
  void Seta2piPDF( const double& a2pi_pdf )           {a2pi_pdf_ = a2pi_pdf;}
  void Setrho0rho0PDF( const double& rho0rho0_pdf )
    {rho0rho0_pdf_ = rho0rho0_pdf;}
  void Setb1piPDF( const double& b1pi_pdf )           {b1pi_pdf_ = b1pi_pdf;}
  void Setrho0pipiPDF( const double& rho0pipi_pdf )
   {rho0pipi_pdf_ = rho0pipi_pdf;}
  void SetpipipipiPDF( const double& pipipipi_pdf )
   {pipipipi_pdf_ = pipipipi_pdf;}

  //Data: 0, MC: 1
  void SetMC( const unsigned int& mc )                {mc_ = mc;}
  //B0B0B event: 0, B+B- event: 1
  void SetBtype( const unsigned int& btype )    {btype_ = btype;}
  void SetTau(double tau)                             {tau_ = tau;}
  void SetDeltam(double delta_m)                      {delta_m_ = delta_m;}

  void SetDeltat( const double& delta_t )             {delta_t_ = delta_t;}

  void SetRecVChisq( const double& recv_chisq )
    {recv_chisq_ = recv_chisq;}
  void SetRecVNdf( const unsigned int& recv_ndf )     {recv_ndf_ = recv_ndf;}
  void SetRecVZerr( const double& recv_zerr )         {recv_zerr_ = recv_zerr;}
  void SetRecVNtrk( const unsigned int& recv_ntrk )   {recv_ntrk_ = recv_ntrk;}

  void SetQ( const int& q )                           {q_ = q;}
  void SetR( const double& r )                        {r_ = r;}
  void SetC( const int& c )                           {c_ = c;}

  void SetTagVChisq( const double& tagv_chisq )
    {tagv_chisq_ = tagv_chisq;}
  void SetTagVNdf( const unsigned int& tagv_ndf )     {tagv_ndf_ = tagv_ndf;}
  void SetTagVZerr( const double& tagv_zerr )         {tagv_zerr_ = tagv_zerr;}
  void SetTagVNtrk( const unsigned int& tagv_ntrk )   {tagv_ntrk_ = tagv_ntrk;}
  void SetTagVIsl( const unsigned int& tagv_isl )     {tagv_isl_ = tagv_isl;}

  void SetAk( const double& a_k )                     {a_k_ = a_k;}
  void SetCk( const double& c_k )                     {c_k_ = c_k;}
  void SetCosTheta( const double& costheta )          {costheta_ = costheta;}

  void SetEfRkRdetRnp( const double& life_pdf )       {life_pdf_ = life_pdf;}
  void SetNormEfRkRdetRnp( const double& int_life_pdf )
    {int_life_pdf_ = int_life_pdf;}
  void SetMfRkRdetRnp( const double& cos_pdf )        {cos_pdf_ = cos_pdf;}
  void SetAfRkRdetRnp( const double& sin_pdf )        {sin_pdf_ = sin_pdf;}

  void SetSigProb( const double& sig_prob )           {sig_prob_ = sig_prob;}
  void SetMisProb( const double& mis_prob )           {mis_prob_ = mis_prob;}
  void SetQQbarProb( const double& qqbar_prob )
    {qqbar_prob_ = qqbar_prob;}
  void SetGB0B0bProb( const double& gb0b0b_prob )
    {gb0b0b_prob_ = gb0b0b_prob;}
  void SetGBpBmProb( const double& gbpbm_prob )
    {gbpbm_prob_ = gbpbm_prob;}
  void SetRB0B0bProb( const double& rb0b0b_prob )
    {rb0b0b_prob_ = rb0b0b_prob;}
  void SetRBpBmProb( const double& rbpbm_prob )
    {rbpbm_prob_ = rbpbm_prob;}
  void Seta2piProb( const double& a2pi_prob )         {a2pi_prob_ = a2pi_prob;}
  void Setrho0rho0Prob( const double& rho0rho0_prob )
    {rho0rho0_prob_ = rho0rho0_prob;}
  void Setb1piProb( const double& b1pi_prob )         {b1pi_prob_ = b1pi_prob;}
  void Setrho0pipiProb( const double& rho0pipi_prob )
    {rho0pipi_prob_ = rho0pipi_prob;}
  void SetpipipipiProb( const double& pipipipi_prob )
    {pipipipi_prob_ = pipipipi_prob;}

 private:
  unsigned int exp_no_;
  unsigned int run_no_;
  unsigned int evt_no_;

  double       mbc_;
  double       benergy_;
  double       de_;
  double       lr_;
  double       fd_;
  double       momega_;
  double       homega_;
  double       mrho0_;

  bool         good_mc_;
  double       momega_gen_;

  double       l0c_;
  double       l2c_;
  double       l0n_;
  double       l2n_;
  double       ctbto_;
  double       ctbz_;
  double       cbz_;
  double       sumpt_;

  double       pi1_ptot_;
  double       pi2_ptot_;
  double       pi3_ptot_;
  double       pi4_ptot_;

  double       pi1_ctheta_;
  double       pi2_ctheta_;
  double       pi3_ctheta_;
  double       pi4_ctheta_;

  double       sig_pdf_;
  double       mis_pdf_;
  double       qqbar_pdf_;
  double       gb0b0b_pdf_;
  double       gbpbm_pdf_;
  double       rb0b0b_pdf_;
  double       rbpbm_pdf_;
  double       a2pi_pdf_;
  double       rho0rho0_pdf_;
  double       b1pi_pdf_;
  double       rho0pipi_pdf_;
  double       pipipipi_pdf_;

  unsigned int mc_;
  unsigned int btype_;
  double       tau_;
  double       delta_m_;

  double       delta_t_;

  double       recv_chisq_;
  unsigned int recv_ndf_;
  double       recv_zerr_;
  unsigned int recv_ntrk_;

  int          q_;
  double       r_;
  int          c_;

  double       tagv_chisq_;
  unsigned int tagv_ndf_;
  double       tagv_zerr_;
  unsigned int tagv_ntrk_;
  unsigned int tagv_isl_;

  double       a_k_;
  double       c_k_;
  double       costheta_;

  double       life_pdf_;
  double       int_life_pdf_;
  double       cos_pdf_;
  double       sin_pdf_;

  double       sig_prob_;
  double       mis_prob_;
  double       qqbar_prob_;
  double       gb0b0b_prob_;
  double       gbpbm_prob_;
  double       rb0b0b_prob_;
  double       rbpbm_prob_;
  double       a2pi_prob_;
  double       rho0rho0_prob_;
  double       b1pi_prob_;
  double       rho0pipi_prob_;
  double       pipipipi_prob_;

  mutable double cache_cos_pdf;
  mutable double cache_sin_pdf;
  mutable double cache_life_pdf;
  mutable double cache_int_life_pdf;
  mutable double cache_lifetime;

 public:
  bool getDTComponents(double lifetime, double& life_pdf, double& int_life_pdf, double& cos_pdf, double& sin_pdf) const {
      if(cache_lifetime != lifetime) return false;
      lifetime = cache_lifetime;
      life_pdf = cache_life_pdf;
      int_life_pdf = cache_int_life_pdf;
      cos_pdf = cache_cos_pdf;
      sin_pdf = cache_sin_pdf;
      return true;
  }
  void setDTComponents(double lifetime, double life_pdf, double int_life_pdf, double cos_pdf, double sin_pdf) const {
      cache_lifetime = lifetime;
      cache_life_pdf = life_pdf;
      cache_int_life_pdf = int_life_pdf;
      cache_cos_pdf = cos_pdf;
      cache_sin_pdf = sin_pdf;
  }

};

void FillPlot( std::ifstream& datafile1, std::ifstream& datafile2,
	       std::vector<Data>& data,
	       const unsigned int& mc_type=0, const unsigned int& B_type=0,
	       const unsigned int& mc_keep=0,
	       const double& ctbto_ul=ctbto_cut,
	       const double& mbc_ll=h_mbc_ll, const double& mbc_ul=h_mbc_ul,
	       const double& de_ll=h_de_ll,  const double& de_ul=h_de_ul,
	       const double& homega_ll=h_homega_ll, const double& homega_ul=h_homega_ul );
void FillPlotGen( std::ifstream& datafile, std::vector<Data>& data );

std::vector<unsigned int> strip_expno( std::vector<Data>& data );
std::vector<double> strip_recv_chisq( std::vector<Data>& data );
std::vector<unsigned int> strip_recv_ndf( std::vector<Data>& data );
std::vector<double> strip_recv_zerr( std::vector<Data>& data );
std::vector<double> strip_tagv_chisq( std::vector<Data>& data );
std::vector<unsigned int> strip_tagv_ndf( std::vector<Data>& data );
std::vector<double> strip_tagv_zerr( std::vector<Data>& data );
std::vector<unsigned int> strip_tagv_isl( std::vector<Data>& data );

std::vector<Data> Add( std::vector<Data>& data1, std::vector<Data>& data2 );

#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

#endif //Data_H_
