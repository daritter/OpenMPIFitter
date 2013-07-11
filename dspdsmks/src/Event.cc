#include "Event.h"
#include <iostream>
#include <iomanip>

size_t dTCache::requests = 0;
size_t dTCache::hits = 0;
size_t dTCache::sets = 0;

void dTCache::print_stats(){
    std::cout << "dT Cache Statistics: "
        << requests << " requests, "
        << hits << " hits, "
        << sets << " sets, ";
    if(requests>0){
        std::ios  state(NULL);
        state.copyfmt(std::cout);
        std::cout << std::fixed << std::setprecision(2) << (100.*hits/requests) << "% efficiency";
        std::cout.copyfmt(state);
    }
    std::cout << std::endl;

    if(hits+sets != requests){
        std::cout << "Warning: hits+sets != requests" << std::endl;
    }
}

void Event::setBranches(TTree* tree, const std::string &bselection){
    std::string prefix("");
#define BADDRESS__(x,var) tree->SetBranchAddress((prefix + x).c_str(),&var)
#define BADDRESS(x) BADDRESS__(#x,x)
#define BODDRESS(g,x) BADDRESS__(#g+"."+#x, g##_##x)
#define BUDDRESS(g,x,y) BADDRESS__(#g+"."+#x, g##_##y)
    BADDRESS(expNo);
    BADDRESS(svdVs);
    BADDRESS(isMC);
    BADDRESS(benergy);
    prefix += bselection + ".";
    BADDRESS(flag);
    BADDRESS(Mbc);
    BADDRESS(dE);
    BADDRESS(m2DspKs);
    BADDRESS(m2DsmKs);
    BADDRESS(cosTheta);
    BADDRESS(deltaZ);
    BODDRESS(tag,ntrk);
    BODDRESS(tag,zerr);
    BODDRESS(tag,chi2);
    BODDRESS(tag,ndf);
    BODDRESS(tag,isL);
    BUDDRESS(tag,flavour,q);
    BUDDRESS(tag,qr,r);
    BODDRESS(vtx,ntrk);
    BODDRESS(vtx,zerr);
    BODDRESS(vtx,chi2);
    BODDRESS(vtx,ndf);
}

void Event::createBranches(TTree* tree, const std::string &bselection){
    std::string prefix("");
#define ADDBRANCH__(name,var,type)   tree->Branch((prefix + name).c_str(),&var,(prefix+name+"/"+#type).c_str())
#define ADDBRANCH(x,type)            ADDBRANCH__(#x,x,type)
#define ADDBRONCH(g,x,type)          ADDBRANCH__(#g+"."+#x,g##_##x,type)
#define ADDBRUNCH(g,x,y,type)        ADDBRANCH__(#g+"."+#x,g##_##y,type)
    ADDBRANCH(expNo,I);
    ADDBRANCH(svdVs,I);
    ADDBRANCH(isMC,O);
    ADDBRANCH(benergy,D);
    prefix += bselection + ".";
    ADDBRANCH(flag,I);
    ADDBRANCH(Mbc,D);
    ADDBRANCH(dE,D);
    ADDBRANCH(m2DspKs,D);
    ADDBRANCH(m2DsmKs,D);
    ADDBRANCH(cosTheta,D);
    ADDBRANCH(deltaZ,D);
    ADDBRONCH(tag,ntrk,I);
    ADDBRONCH(tag,zerr,D);
    ADDBRONCH(tag,chi2,D);
    ADDBRONCH(tag,ndf,D);
    ADDBRONCH(tag,isL,D);
    ADDBRUNCH(tag,flavour,q,I);
    ADDBRUNCH(tag,qr,r,D);
    ADDBRONCH(vtx,ntrk,I);
    ADDBRONCH(vtx,zerr,D);
    ADDBRONCH(vtx,chi2,D);
    ADDBRONCH(vtx,ndf,I);
}
