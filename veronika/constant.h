#ifndef Constant_H_
#define Constant_H_

// Particle masses
const double m_pi = 0.13957018;
const double m_rho = 0.77549;
const double m_a2pi_mc = 1.3183;
const double m_b1pi_mc = 1.2295;
const double m_omega_mc = 0.782;

// Branching fractions
const double bf_rho0rho0 = 0.73e-06;
const double bf_rho0rho0_err = 0.28e-06;
const double bf_b1pi = 10.9e-06*0.0153;
const double bf_b1pi_err =
  bf_b1pi*
  sqrt(((1.5e-06/10.9e-06)*(1.5e-06/10.9e-06)) +
       ((0.0013/0.0153)*(0.0013/0.0153)));
const double bf_rho0pipi = 4.35e-06;
const double bf_pipipipi = 9.65e-06;

// Selection criteria
const double ctbto_cut = 0.9;

// Fit region
const double h_mbc_ll = 5.25;
const double h_mbc_ul = 5.30;
const double h_de_ll = -0.15;
const double h_de_ul = 0.1;
const double h_fd_ll = -10.;
const double h_fd_ul = 10.;
const double h_momega_ll = 0.73;
const double h_momega_ul = 0.83;
const double h_homega_ul = 1.0;
const double h_homega_ll = -1.0;
const double h_lr_ll = 0.2;
const double h_lr_ul = 1.;

//Binning
const unsigned int bins_de  = 100;
const unsigned int bins_fd  = 25;
const unsigned int bins_momega = 100;
const unsigned int bins_homega = 50;
const unsigned int bins_dt = 50;
const unsigned int bins_mbc = 100;

// Signal region
const double sb_de_ll = -0.04;
const double sb_de_ul = -sb_de_ll;
const double sb_fd_ll = 0.0;
const double sb_dt_ll = -7.5;
const double sb_dt_ul = -sb_dt_ll;

// Sideband region
const double side_mbc_ll = 5.20;
const double side_mbc_ul = 5.25;
const double side_de_ll = 0.05;
const double side_de_ul = 0.2;

// MC scale factors
const double sf_generic = 1.0/6.0;
const double sf_rareb = 1.0/50.0;

// NBB
const double NBB_svd1 = 151.961e+06;
const double NBB_err_svd1 = 1.241e+06;
const double NBB_svd2 = 619.620e+06;
const double NBB_err_svd2 = 9.441e+06;
const double NBB_mc_pb = 500000.0;

// Efficiency
const double omegaKs_eff_svd1 = 109724.0/1000000.;
const double omegaKs_eff_err_svd1 = sqrt(109724.0)/1000000.;
const double omegaKs_eff_svd2 = 141283.0/1000000.;
const double omegaKs_eff_err_svd2 = sqrt(141283.0)/1000000.;
const double eff_w = 0.89401;
const double misrec_frac_svd1 = 1331./(109724.);
const double misrec_frac_svd2 = 1703./(141283.);

//TCPV constants
const double b0_lifetime = 1.525;
const double delta_m = 0.507;
const double b0_lifetime_mc = 1.53439;
const double delta_m_mc = 0.507;

// TCPV cuts
const double dt_cut = 70.0;
const double h_cut = 50.0;
const double sigmaz_mult_cut = 0.02;
const double sigmaz_sngl_cut = 0.05;

const double h_dt_ll = -dt_cut;
const double h_dt_ul = dt_cut;

const double h_q_ll = -1.;
const double h_q_ul = 1.;

const double tcpv_dt_ll = -7.5;
const double tcpv_dt_ul = 7.5;

#endif //Constant_H_
