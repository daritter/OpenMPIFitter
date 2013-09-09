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
        BinEfficiency(int index):index(index),channel_count{{{0}}},count{0} {
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

            slowPi[0].init(TRK_SVD1, true);
            slowPi[1].init(TRK_SVD2, true);
            slowPi0[0].init(PI0_SVD1, false);
            slowPi0[1].init(PI0_SVD2, false);
            fastPi0[0].init(PI0_SVD1, false);
            fastPi0[1].init(PI0_SVD2, false);
            fastTracks[0].init(TRK_SVD1, true);
            fastTracks[1].init(TRK_SVD2, true);
            //init(slowPi,   "slowPi",  false, 0.0);
        }

        ~BinEfficiency() {}

        void init(Belle::KID_eff_06* kid_eff, const std::string name, bool K, double pid=0.1){
            kid_eff[0].init(pid, K?1:3, (prefix + name + "_svd1").c_str(), SVD1_TABLE);
            kid_eff[1].init(pid, K?1:3, (prefix + name + "_svd2").c_str(), SVD2_TABLE);
        }

        int addDs(int svd, int channel, int* pdg, double* plab, double* costheta){
            //Use the D channel numbers, not D*
            if(channel>5) channel -= 4;
            channel -= 2;
            /*std::cout << channel << ": ";
            for(int i=0; i<5; ++i){
                std::cout << "(" << pdg[i] << "," << plab[i] << "," << costheta[i] << "), ";
            }
            std::cout << std::endl;*/
            switch(channel){
                case 0:
                    assert(fabs(pdg[0]) == 321);
                    assert(fabs(pdg[1]) == 211);
                    assert(fabs(pdg[2]) == 0);
                    assert(fabs(pdg[3]) == 0);
                    assert(fabs(pdg[4]) == 211);
                    DKp_K[svd].addtrack(plab[0], costheta[0]);
                    DKp_p[svd].addtrack(plab[1], costheta[1]);
                    fastTracks[svd].add(plab[0]);
                    fastTracks[svd].add(plab[1]);
                    break;
                case 1:
                    assert(fabs(pdg[0]) == 321);
                    assert(fabs(pdg[1]) == 211);
                    assert(fabs(pdg[2]) == 111);
                    assert(fabs(pdg[3]) == 0);
                    assert(fabs(pdg[4]) == 211);
                    DKpp0_K[svd].addtrack(plab[0], costheta[0]);
                    DKpp0_p[svd].addtrack(plab[1], costheta[1]);
                    fastPi0[svd].add(plab[2]);
                    fastTracks[svd].add(plab[0]);
                    fastTracks[svd].add(plab[1]);
                    break;
                case 2:
                    assert(fabs(pdg[0]) == 211);
                    assert(fabs(pdg[1]) == 211);
                    assert(fabs(pdg[2]) == 321);
                    assert(fabs(pdg[3]) == 211);
                    assert(fabs(pdg[4]) == 211);
                    DKppp_p1[svd].addtrack(plab[0], costheta[0]);
                    DKppp_p2[svd].addtrack(plab[1], costheta[1]);
                    DKppp_K[svd].addtrack(plab[2], costheta[2]);
                    DKppp_p3[svd].addtrack(plab[3], costheta[3]);
                    fastTracks[svd].add(plab[0]);
                    fastTracks[svd].add(plab[1]);
                    fastTracks[svd].add(plab[2]);
                    fastTracks[svd].add(plab[3]);
                    break;
                case 3:
                    assert(fabs(pdg[0]) == 211);
                    assert(fabs(pdg[1]) == 211);
                    assert(fabs(pdg[2]) == 321);
                    assert(fabs(pdg[3]) == 0);
                    assert(fabs(pdg[4]) == 111);
                    DKpp_p1[svd].addtrack(plab[0], costheta[0]);
                    DKpp_p2[svd].addtrack(plab[1], costheta[1]);
                    DKpp_K[svd].addtrack(plab[2], costheta[2]);
                    fastTracks[svd].add(plab[0]);
                    fastTracks[svd].add(plab[1]);
                    fastTracks[svd].add(plab[2]);
                    break;
                default:
                    std::cerr << "Unknown channel: " << channel << std::endl;
                    std::abort();
            }
            if(fabs(pdg[4]) == 211) {
                slowPi[svd].add(plab[4]);
            }else{
                slowPi0[svd].add(plab[4]);
            }
            return channel;
        }

        void add(int svd, int channelP, int channelM, int* pdg, double* plab, double* costheta) {
            /*std::cout << channelP << ", " << channelM << ": ";
            for(int i=0; i<11; ++i){
                std::cout << pdg[i] << ", ";
            }
            std::cout << std::endl;*/
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

        void calc_eff_err(double& eff, double& err, std::vector<Belle::KID_eff_06*> kid1, const std::vector<Belle::KID_eff_06*> &kid2, const std::vector<Belle::KID_eff_06*> &kid3, bool usemc=false){
            kid1.insert(end(kid1), begin(kid2), end(kid2));
            kid1.insert(end(kid1), begin(kid3), end(kid3));
            eff = 1.0;
            err = 0;
            for(size_t i=0; i<kid1.size(); ++i){
                const auto& kid = kid1[i];
                eff *= usemc?kid->ratio_ref():kid->ratio();
                double tmp(1.);
                for(size_t j=0; j<kid1.size(); ++j){
                    if(i==j) continue;
                    const auto& kid2 = kid1[j];
                    tmp *= usemc?kid->ratio_ref():kid2->ratio();
                }
                err += tmp * kid->ratio_error();
            }
        }


        void save(std::ostream &out) {
            std::vector<std::string> channel_names{
                "D⁰ → K⁻π⁺    ",
                "D⁰ → K⁻π⁺π⁰  ",
                "D⁰ → K⁻π⁺π⁺π⁻",
                "D⁺ → K⁻π⁺π⁺  "
            };
            boost::format value_error("%.4f ± %.4f (%.2f%%) -- weight = %.4f");
            out << "Efficiencies for Bin " << index << ":" << std::endl;
            out.precision(5);
            for(int svd=0; svd<2; ++svd){
                out << " =============================== SVD " << (svd+1)
                    << " =============================== " << std::endl;
                std::vector<Belle::KID_eff_06*> channels[4] = {
                    {&DKp_K[svd], &DKp_p[svd]},
                    {&DKpp0_K[svd], &DKpp0_p[svd]},
                    {&DKppp_K[svd], &DKppp_p1[svd], &DKppp_p2[svd], &DKppp_p3[svd]},
                    {&DKpp_K[svd], &DKpp_p1[svd], &DKpp_p2[svd]}
                };
                std::vector<Belle::KID_eff_06*> Km;
                if(!ks_eff[svd].get_totalNum()>0){
                    Km.push_back(&km_eff[svd]);
                }

                for(int i=0; i<4; ++i){
                    std::cout << "⇒ KID_eff values for " << channel_names[i] << std::endl;
                    for(auto* kid: channels[i]){
                        kid->calculate();
                        kid->dump();
                    }
                }
                double eff[4][4] = {{0}};
                double err[4][4] = {{0}};
                double total_eff = 0;
                double total_err = 0;
                double total_weight = 0;
                for(int i=0; i<4; ++i){
                    for(int j=0; j<4; ++j){
                        double weight = 1.0 * channel_count[svd][i][j] / count[svd];
                        calc_eff_err(eff[i][j], err[i][j], channels[i], channels[j], Km);
                        out << "Channel " << channel_names[i] << ", " << channel_names[j] << ": "
                            << value_error % eff[i][j] %  err[i][j] % (err[i][j]/eff[i][j]*100) % weight
                            << " (" << channel_count[svd][i][j] << "/" << count[svd] << ")" << std::endl;
                        total_eff += weight*eff[i][j];
                        total_err += weight*weight*err[i][j]*err[i][j];
                        total_weight += weight;
                    }
                }
                total_err = std::sqrt(total_err);
                out << "≣≣≣≣ KID Efficiency correction: "
                    << value_error % total_eff % total_err % (total_err/total_eff*100) % total_weight
                    << std::endl;
                if(ks_eff[svd].get_totalNum()>0){
                    double ks_err = ks_eff[svd].total_errFactor() /  ks_eff[svd].total_effFactor();
                    ks_err = std::sqrt((ks_err)*(ks_err) + 0.006*0.006);
                    ks_err *= ks_eff[svd].total_effFactor();
                    out << "≣≣≣≣  Ks Efficiency correction: "
                        << value_error % ks_eff[svd].total_effFactor()
                        % ks_err
                        % (ks_err/ks_eff[svd].total_effFactor()*100) % 1.0
                       << std::endl;
                }
                slowPi[svd].calculate();
                out << "≣≣≣≣ slowPi Efficiency correction: "
                    << value_error % slowPi[svd].get_eff() % slowPi[svd].get_err()
                    % (slowPi[svd].get_err()/slowPi[svd].get_eff()*100) % 1.0
                    << std::endl;
                slowPi0[svd].calculate();
                out << "≣≣≣≣ slowPi0 Efficiency correction: "
                    << value_error % slowPi0[svd].get_eff() % slowPi0[svd].get_err()
                    % (slowPi0[svd].get_err()/slowPi0[svd].get_eff()*100) % 1.0
                    << std::endl;
                fastPi0[svd].calculate();
                out << "≣≣≣≣ fastPi0 Efficiency correction: "
                    << value_error % fastPi0[svd].get_eff() % fastPi0[svd].get_err()
                    % (fastPi0[svd].get_err()/fastPi0[svd].get_eff()*100) % 1.0
                    << std::endl;
                fastTracks[svd].calculate();
                out << "≣≣≣≣ fastTracks Efficiency correction: "
                    << value_error % fastTracks[svd].get_eff() % fastTracks[svd].get_err()
                    % (fastTracks[svd].get_err()/fastTracks[svd].get_eff()*100) % 1.0
                    << std::endl;
                out << std::endl;
            }
        }


    protected:
        int index;
        std::string prefix;
        unsigned int channel_count[2][4][4];
        unsigned int count[2];

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
        TrackEfficiency slowPi[2];
        TrackEfficiency slowPi0[2];

        TrackEfficiency fastTracks[2];
        TrackEfficiency fastPi0[2];
};

#endif
