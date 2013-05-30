#include <Parameters.h>
#include "dspdsmks/DspDsmKs.h"
#include "gtest/gtest.h"
#include "TChain.h"
#include <iostream>
#include <fstream>

namespace {
    inline double getStep(int i, int N, double min, double max){
        return 1.0*i/(N-1) * (max-min) + min;
    }
}

#define step(var,min,max,N) for(double i##var=0, var=min; i##var<N; i##var+=1, var=1.0*i##var/(N-1)*(max-min)+min)

TEST(DeltaT,NormWithCache){
    const double dtmin = -70;
    const double dtmax =  70;
    const int dtsteps = 10001;
    Range range_mBC("Mbc", 5.24,5.3);
    Range range_dE("dE", -0.15,0.1);
    DeltaTPDF dtpdf(Range("dT",-70,70), 0);
    std::vector<double> par({1.53, 0, 0, 0});
    dtpdf.setParameters(0,1,2,3,-1,-1, 0);

    TChain* chain = new TChain("B0");
    chain->AddFile("/home/iwsatlas1/ritter/belle/DspDsmKs/skim/ddk-on_resonance.root");
    Event event;
    event.setBranches(chain, "bestLHsig");
    for(int svd=0; svd<2; ++svd){
        int i(0);
        do{
            chain->GetEntry(i++);
        }while(event.svdVs != svd || !event.calculateValues(true) || !range_mBC(event.Mbc) || !range_dE(event.dE));

        step(Jc,-1,1,5){ step(Js1,-1,1,5) { step(Js2,-1,1,5){
            par[1] = Jc; par[2] = Js1; par[3] = Js2;
            for(event.eta=-1; event.eta<2; event.eta+=2){
                long double norm(0);
                step(dt,dtmin,dtmax,dtsteps){
                    event.deltaT = dt;
                    for(event.tag_q=-1; event.tag_q<2; event.tag_q+=2){
                        norm += dtpdf(event, par);
                    }
                }
                norm *= (dtmax-dtmin)/(dtsteps-1);
                EXPECT_FLOAT_EQ(norm, 0.5) << "Jc: " << Jc << ", Js1: " << Js1 << ", Js2: " << Js2;
            }
        }}}
        std::cout << "svd " << (svd+1) << " done" << std::endl;
    }
    dTCache::print_stats();
    delete chain;
}

TEST(BkgT,NormWithCache){
    const double dtmin = -70;
    const double dtmax =  70;
    const int dtsteps = 10001;
    Range range_mBC("Mbc", 5.24,5.3);
    Range range_dE("dE", -0.15,0.1);

    BkgTPDF dtpdf(Range("dT",-70,70), 0);
    BBarPDF::init_deltaT(dtpdf,false);

    Parameters params;
    {
        std::ifstream file("ctrl-out.par");
        file >> params;
        file.close();
    }
    std::vector<double> par = params.getValues();

    TChain* chain = new TChain("B0");
    chain->AddFile("/home/iwsatlas1/ritter/belle/DspDsmKs/skim/ddk-on_resonance.root");
    Event event;
    event.setBranches(chain, "bestLHsig");
    for(int svd=0; svd<2; ++svd){
        int i(0);
        do{
            chain->GetEntry(i++);
        }while(event.svdVs != svd || !event.calculateValues(true) || !range_mBC(event.Mbc) || !range_dE(event.dE));

        for(unsigned int p=PAR::bkg_dt_blifetime; p<PAR::bkg_svd2_dt_acp6; ++p){
            double v = params[p].value();
            double e = params[p].error();
            step(pvalue, v-e, v+e, 10){
                par[p] = pvalue;
                long double norm(0);
                step(dt,dtmin,dtmax,dtsteps){
                    event.deltaT = dt;
                    for(event.tag_q=-1; event.tag_q<2; event.tag_q+=2){
                        for(event.eta=-1; event.eta<2; event.eta+=2){
                            norm += dtpdf(event, par);
                        }
                    }
                }
                norm *= (dtmax-dtmin)/(dtsteps-1);
                ASSERT_FLOAT_EQ(norm, 1.) << "svd: " << svd << ", " << params[p].name() << ": " << pvalue;
            }
            std::cout << "svd" << (svd+1) << ", " <<  params[p].name() << " done" << std::endl;
            par[p] = v;
        }
    }
    delete chain;
    dTCache::print_stats();
}
