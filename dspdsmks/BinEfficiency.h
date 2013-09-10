#ifndef MPIFitter_BinEfficiency_h
#define MPIFitter_BinEfficiency_h

#include "Ks_eff.h"
namespace Belle {
#include "kid_eff_06s.h"
}
#include "TrackEfficiency.h"
#include <boost/format.hpp>
#include <cassert>
#include <cstdlib>

#define SVD1_TABLE "dspdsmks/kid/kideff-2006-svd1-all.dat"
#define SVD2_TABLE "dspdsmks/kid/kideff-2010.dat"
#define TRK_SVD1 "dspdsmks/trk/slowPi-SVD1.txt"
#define TRK_SVD2 "dspdsmks/trk/slowPi-SVD2.txt"
#define PI0_SVD1 "dspdsmks/trk/slowPi0-SVD1.txt"
#define PI0_SVD2 "dspdsmks/trk/slowPi0-SVD2.txt"

class BinEfficiency {
    public:
        BinEfficiency(int index):index(index),channel_count{{{0}}},count{0},value_error("%.4f ± %.4f (%5.2f%%)") {
            prefix = (boost::format("bin%1%_") % index).str();

            init(km_eff,   "Km",      true);
            init(DKp_K,    "DKp_K",   true);
            init(DKp_p,    "DKp_p",   false);
            init(DKpp0_K,  "DKpp0_K", true);
            init(DKpp0_p,  "DKpp0_p", false);
            init(DKppp_K,  "DKpp_K",  true);
            init(DKppp_p1, "DKpp_p1", false);
            init(DKppp_p2, "DKpp_p2", false);
            init(DKppp_p3, "DKpp_p3", false);
            init(DKpp_K,   "DKpp_K",  true);
            init(DKpp_p1,  "DKpp_p1", false);
            init(DKpp_p2,  "DKpp_p2", false);

            init(trk_slowPi, true);
            init(trk_slowPi0, false);
            init(trk_Km_eff, true);
            init(trk_DKp_K, true);
            init(trk_DKp_p, true);
            init(trk_DKpp0_K, true);
            init(trk_DKpp0_p, true);
            init(trk_DKpp0_p0, false);
            init(trk_DKppp_K, true);
            init(trk_DKppp_p1, true);
            init(trk_DKppp_p2, true);
            init(trk_DKppp_p3, true);
            init(trk_DKpp_K, true);
            init(trk_DKpp_p1, true);
            init(trk_DKpp_p2, true);
        }

        ~BinEfficiency() {}

        void init(Belle::KID_eff_06* kid_eff, const std::string name, bool K, double pid=0.1){
            kid_eff[0].init(pid, K?1:3, (prefix + name + "_svd1").c_str(), SVD1_TABLE);
            kid_eff[1].init(pid, K?1:3, (prefix + name + "_svd2").c_str(), SVD2_TABLE);
        }

        void init(TrackEfficiency* trk_eff, bool charged){
            trk_eff[0].init(charged?TRK_SVD1:PI0_SVD1, charged);
            trk_eff[1].init(charged?TRK_SVD2:PI0_SVD2, charged);
        }

        int addDs(int svd, int channel, int* pdg, double* plab, double* costheta){
            //Use the D channel numbers, not D*
            if(channel>5) channel -= 4;
            channel -= 2;
            switch(channel){
                case 0:
                    assert(fabs(pdg[0]) == 321);
                    assert(fabs(pdg[1]) == 211);
                    assert(fabs(pdg[2]) == 0);
                    assert(fabs(pdg[3]) == 0);
                    DKp_K[svd].addtrack(plab[0], costheta[0]);
                    DKp_p[svd].addtrack(plab[1], costheta[1]);
                    trk_DKp_K[svd].add(plab[0]);
                    trk_DKp_p[svd].add(plab[1]);
                    break;
                case 1:
                    assert(fabs(pdg[0]) == 321);
                    assert(fabs(pdg[1]) == 211);
                    assert(fabs(pdg[2]) == 111);
                    assert(fabs(pdg[3]) == 0);
                    DKpp0_K[svd].addtrack(plab[0], costheta[0]);
                    DKpp0_p[svd].addtrack(plab[1], costheta[1]);
                    trk_DKpp0_K[svd].add(plab[0]);
                    trk_DKpp0_p[svd].add(plab[1]);
                    trk_DKpp0_p0[svd].add(plab[2]);
                    break;
                case 2:
                    assert(fabs(pdg[0]) == 211);
                    assert(fabs(pdg[1]) == 211);
                    assert(fabs(pdg[2]) == 321);
                    assert(fabs(pdg[3]) == 211);
                    DKppp_p1[svd].addtrack(plab[0], costheta[0]);
                    DKppp_p2[svd].addtrack(plab[1], costheta[1]);
                    DKppp_K[svd].addtrack(plab[2], costheta[2]);
                    DKppp_p3[svd].addtrack(plab[3], costheta[3]);
                    trk_DKppp_p1[svd].add(plab[0]);
                    trk_DKppp_p2[svd].add(plab[1]);
                    trk_DKppp_K[svd].add(plab[2]);
                    trk_DKppp_p3[svd].add(plab[3]);
                    break;
                case 3:
                    assert(fabs(pdg[0]) == 211);
                    assert(fabs(pdg[1]) == 211);
                    assert(fabs(pdg[2]) == 321);
                    assert(fabs(pdg[3]) == 0);
                    DKpp_p1[svd].addtrack(plab[0], costheta[0]);
                    DKpp_p2[svd].addtrack(plab[1], costheta[1]);
                    DKpp_K[svd].addtrack(plab[2], costheta[2]);
                    trk_DKpp_p1[svd].add(plab[0]);
                    trk_DKpp_p2[svd].add(plab[1]);
                    trk_DKpp_K[svd].add(plab[2]);
                    break;
                default:
                    std::cerr << "Unknown channel: " << channel << std::endl;
                    std::abort();
            }
            if(fabs(pdg[4]) == 211) {
                trk_slowPi[svd].add(plab[4]);
            }else if(fabs(pdg[4]) == 111){
                trk_slowPi0[svd].add(plab[4]);
            }
            return channel;
        }

        void add(int svd, int channelP, int channelM, int* pdg, double* plab, double* costheta) {
            int c1 = addDs(svd, channelP, pdg, plab, costheta);
            int c2 = addDs(svd, channelM, pdg + 5, plab + 5, costheta + 5);
            ++channel_count[svd][c1][c2];
            ++count[svd];

            //Ks_eff:
            if(fabs(pdg[10]) == 310 || fabs(pdg[10]) == 311){
                ks_eff[svd].add(plab[10], costheta[10]);
            }else{
                km_eff[svd].addtrack(plab[10], costheta[10]);
            }
        }

        template<class T> void calc_eff_err(double& eff, double& err, std::vector<T*> v1, const std::vector<T*> &v2, const std::vector<T*> &v3, bool use_ref=false, bool conservative=false){
            // Calculate the efficiency and error for one channel. We have two
            // d(*)s and possibly a K- (in the control channel) so lets add up
            // all the tracks into one vector
            v1.insert(end(v1), begin(v2), end(v2));
            v1.insert(end(v1), begin(v3), end(v3));

            // No let's calculate the combined efficiency and error
            // eta = prod_i eta_i
            // sigma = sum_i (prod_j^(j!=i) eta_j) sigma_i
            eff = 1.0;
            err = 0;
            for(size_t i=0; i<v1.size(); ++i){
                const auto& e1 = v1[i];
                eff *= use_ref?e1->ratio_ref():e1->ratio();
                double tmp(1.);
                for(size_t j=0; j<v1.size(); ++j){
                    if(i==j) continue;
                    const auto& e2 = v1[j];
                    tmp *= use_ref?e2->ratio_ref():e2->ratio();
                }
                double e_err = e1->ratio_error();
                if(conservative){
                    const double d = e1->ratio() - e1->ratio_ref();
                    e_err = std::sqrt(e_err*e_err + d*d);
                }
                err += tmp * e_err;
            }
        }

        template<class T> void calculate_total(std::ostream &out, int svd, double &eff, double &err, const std::vector<T*>* d1, const std::vector<T*>* d2, const std::vector<T*> &k, const std::string &name){
            static std::vector<std::string> channel_names{
                "D⁰ → K⁻π⁺    ",
                "D⁰ → K⁻π⁺π⁰  ",
                "D⁰ → K⁻π⁺π⁺π⁻",
                "D⁺ → K⁻π⁺π⁺  "
            };
            //Calculate efficency and error for all channels and make a
            //weighted average of the efficiency and the relative errors
            for(int i=0; i<4; ++i){
                for(int j=0; j<4; ++j){
                    double weight = 1.0 * channel_count[svd][i][j] / count[svd];
                    if(weight == 0) continue;

                    double ch_eff(0), ch_err(0);
                    calc_eff_err(ch_eff, ch_err, d1[i], d2[j], k);
                    if(ch_err>0){
                        out << name + " " + channel_names[i] + ", " + channel_names[j] << ": "
                            << value_error % ch_eff %  ch_err % (ch_err/ch_eff*100)
                            << boost::format(" -- weight = %.4f") % weight << std::endl;
                    }
                    eff += weight*ch_eff;
                    err += weight*ch_err/ch_eff;
                }
            }
            out << std::endl;
            err *= eff;
        }

        void overall(double &eff, double &err, const std::vector<double> &efficiencies, const std::vector<double> &errors){
            eff = 1;
            err = 0;
            for(size_t i=0; i<efficiencies.size(); ++i){
                eff *= efficiencies[i];
                const double rel = errors[i]/efficiencies[i];
                err += rel*rel;
            }
            err = sqrt(err)*eff;
        }

        void save(std::ostream &out) {
            auto print = [&](double tot_eff, double tot_err, const std::string &name){
                out << "≣≣≣≣ " << name << " efficiency correction: "
                    << value_error % tot_eff % tot_err % (tot_err/tot_eff*100)
                    << std::endl;
            };
            out << "Efficiencies for Bin " << index << ":" << std::endl;
            for(int svd=0; svd<2; ++svd){
                out << "=================================== SVD " << (svd+1)
                    << " ===================================" << std::endl;

                //kid efficienies used in the different decay channels
                std::vector<Belle::KID_eff_06*> kid_channels[4] = {
                    {&DKp_K[svd], &DKp_p[svd]},
                    {&DKpp0_K[svd], &DKpp0_p[svd]},
                    {&DKppp_K[svd], &DKppp_p1[svd], &DKppp_p2[svd], &DKppp_p3[svd]},
                    {&DKpp_K[svd], &DKpp_p1[svd], &DKpp_p2[svd]}
                };
                std::vector<Belle::KID_eff_06*> kid_Km;
                //trk efficiencies used in the different decay channels
                std::vector<TrackEfficiency*> trk_channels[4] = {
                    {&trk_DKp_K[svd], &trk_DKp_p[svd], &trk_slowPi[svd]},
                    {&trk_DKpp0_K[svd], &trk_DKpp0_p[svd], &trk_slowPi[svd]},
                    {&trk_DKppp_K[svd], &trk_DKppp_p1[svd], &trk_DKppp_p2[svd], &trk_DKppp_p3[svd], &trk_slowPi[svd]},
                    {&trk_DKpp_K[svd], &trk_DKpp_p1[svd], &trk_DKpp_p2[svd]}
                };
                //for ctrl channel the second is a D and does not have slow pi
                std::vector<TrackEfficiency*> trk_channels_ctrl[4] = {
                    {&trk_DKp_K[svd], &trk_DKp_p[svd]},
                    {&trk_DKpp0_K[svd], &trk_DKpp0_p[svd]},
                    {&trk_DKppp_K[svd], &trk_DKppp_p1[svd], &trk_DKppp_p2[svd], &trk_DKppp_p3[svd]},
                    {&trk_DKpp_K[svd], &trk_DKpp_p1[svd], &trk_DKpp_p2[svd]}
                };
                std::vector<TrackEfficiency*> trk_Km;
                //pi0 efficiencies used in the different decay channels
                std::vector<TrackEfficiency*> pi0_channels[4] = {
                    {},
                    {&trk_DKpp0_p0[svd]},
                    {},
                    {&trk_slowPi0[svd]}
                };
                std::vector<TrackEfficiency*> pi0_Km;

                bool ctrl=false;
                if(!ks_eff[svd].get_totalNum()>0){
                    kid_Km.push_back(&km_eff[svd]);
                    trk_Km.push_back(&trk_Km_eff[svd]);
                    ctrl = true;
                }

                //Calculate all used efficiencies
                for(int i=0; i<4; ++i){
                    //std::cout << "⇒ KID_eff values for " << channel_names[i] << std::endl;
                    for(auto* kid: kid_channels[i]){ kid->calculate(); }
                    for(auto* trk: trk_channels[i]){ trk->calculate(); }
                    for(auto* trk: pi0_channels[i]){ trk->calculate(); }
                }
                for(auto* trk: trk_Km){ trk->calculate(); }
                for(auto* kid: kid_Km){ kid->calculate(); }

                double kid_total_eff(0), kid_total_err(0);
                double trk_total_eff(0), trk_total_err(0);
                double pi0_total_eff(0), pi0_total_err(0);
                double k0s_total_eff(1), k0s_total_err(0);
                double complete_eff(1), complete_err(0);

                calculate_total(out, svd, kid_total_eff, kid_total_err, kid_channels, kid_channels, kid_Km, "kid");
                calculate_total(out, svd, trk_total_eff, trk_total_err, trk_channels, ctrl?trk_channels_ctrl:trk_channels, trk_Km, "trk");
                calculate_total(out, svd, pi0_total_eff, pi0_total_err, pi0_channels, pi0_channels, pi0_Km, "pi0");

                print(kid_total_eff, kid_total_err, "kid");
                print(trk_total_eff, trk_total_err, "trk");
                print(pi0_total_eff, pi0_total_err, "pi0");

                if(ks_eff[svd].get_totalNum()>0){
                    k0s_total_eff = ks_eff[svd].total_effFactor();
                    k0s_total_err = ks_eff[svd].total_errFactor() /  ks_eff[svd].total_effFactor();
                    k0s_total_err = std::sqrt((k0s_total_err)*(k0s_total_err) + 0.006*0.006);
                    k0s_total_err *= k0s_total_eff;
                    print(k0s_total_eff, k0s_total_err, "K0s");
                }

                //Calculate the overall efficiency&error (adding relative error quadratically)
                overall(complete_eff, complete_err,
                        {kid_total_eff, trk_total_eff, pi0_total_eff, k0s_total_eff},
                        {kid_total_err, trk_total_err, pi0_total_err, k0s_total_err}
                );

                print(complete_eff, complete_err, "all");

                out << std::endl;
            }
        }


    protected:
        int index;
        std::string prefix;
        unsigned int channel_count[2][4][4];
        unsigned int count[2];
        boost::format value_error;

        Belle::Ks_eff ks_eff[2];
        Belle::KID_eff_06 km_eff[2];

        //Channel 2: D -> KmPp
        Belle::KID_eff_06 DKp_K[2];
        Belle::KID_eff_06 DKp_p[2];

        //Channel 3: D -> KmPpP0
        Belle::KID_eff_06 DKpp0_K[2];
        Belle::KID_eff_06 DKpp0_p[2];

        //Channel 4: D -> KmPpPpPm
        Belle::KID_eff_06 DKppp_K[2];
        Belle::KID_eff_06 DKppp_p1[2];
        Belle::KID_eff_06 DKppp_p2[2];
        Belle::KID_eff_06 DKppp_p3[2];

        //Channel 5: D -> KmPpPp
        Belle::KID_eff_06 DKpp_K[2];
        Belle::KID_eff_06 DKpp_p1[2];
        Belle::KID_eff_06 DKpp_p2[2];

        //Slow pions for all channels the same
        //KID_eff_06 slowPi[2];
        TrackEfficiency trk_slowPi[2];
        TrackEfficiency trk_slowPi0[2];

        TrackEfficiency trk_Km_eff[2];

        //Channel 2: D -> KmPp
        TrackEfficiency trk_DKp_K[2];
        TrackEfficiency trk_DKp_p[2];

        //Channel 3: D -> KmPpP0
        TrackEfficiency trk_DKpp0_K[2];
        TrackEfficiency trk_DKpp0_p[2];
        TrackEfficiency trk_DKpp0_p0[2];

        //Channel 4: D -> KmPpPpPm
        TrackEfficiency trk_DKppp_K[2];
        TrackEfficiency trk_DKppp_p1[2];
        TrackEfficiency trk_DKppp_p2[2];
        TrackEfficiency trk_DKppp_p3[2];

        //Channel 5: D -> KmPpPp
        TrackEfficiency trk_DKpp_K[2];
        TrackEfficiency trk_DKpp_p1[2];
        TrackEfficiency trk_DKpp_p2[2];
};

#endif
