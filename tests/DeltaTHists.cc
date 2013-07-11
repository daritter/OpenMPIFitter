#include "gtest/gtest.h"
#include "dspdsmks/DeltaTHists.h"
#include "TFile.h"
#include "TF1.h"
#include <iostream>

TEST(DeltaTHists, Fill){
    //TFile *f = new TFile("testhists.root","RECREATE");
    DeltaTHists dthists(10,-0.5,9.5, "test");
    Event e;
    for(e.svdVs=0; e.svdVs<2; ++e.svdVs){
        for(e.rbin=0; e.rbin<7; ++e.rbin){
            dthists.fill(e, [](const Event &e){
                    return e.deltaT * (DeltaTHists::NHISTS_PER_SVD*e.svdVs + DeltaTHists::NHISTS_PER_RBIN*e.rbin + (e.tag_q+1) + (e.eta+1)/2);
                    }, true);
        }
    }

    DeltaTHists dthists2(10,-0.5,9.5, "test2");
    dthists2.recieve(dthists.send());
    for(int svd=0; svd<2; ++svd){
        for(int rbin=0; rbin<7; ++rbin){
            ASSERT_EQ(dthists2.yield(svd, rbin), dthists.yield(svd, rbin));
        }
    }
    for(int i=0; i<DeltaTHists::NHISTS; ++i){
        TH1D* h1 = dthists[i];
        TH1D* h2 = dthists2[i];
        for(int bin=0; bin<h1->GetSize(); ++bin){
            ASSERT_EQ(h1->GetBinContent(bin), h2->GetBinContent(bin));
        }
    }
    //f->Write();
    //f->Close();
}

TEST(DeltaTHists, Norm){
    //TFile *f = new TFile("testhists.root","RECREATE");
    DeltaTHists dthists(5000,-10,10, "gauss");
    Event e;
    e.svdVs = 0;
    e.rbin = 0;
    TF1 gauss("g","gausn",-10,10);
    gauss.SetParameters(1,0,1);
    for(int i=0; i<100; ++i){
        dthists.fill(e, [&](const Event &e){ return gauss(e.deltaT - i*5./100.); }, true);
    }
    DeltaTHists dthists2(5000,-10,10);
    dthists2.recieve(dthists.send());
    dthists.finalize(true);
    dthists2.finalize(true,[](int,int){return 2.;});
    for(int i=0; i<DeltaTHists::NHISTS; ++i){
        ASSERT_FLOAT_EQ(dthists[i]->Integral(), (i<4)?1:0);
        ASSERT_FLOAT_EQ(dthists2[i]->Integral(), (i<4)?2:0);
    }
    //f->Write();
    //f->Close();
}
