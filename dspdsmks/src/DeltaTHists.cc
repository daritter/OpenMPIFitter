#include "DeltaTHists.h"
#include <boost/format.hpp>

namespace {
    /**
     * Set a histogram from a std::vector, starting at the given index and
     * increase the index by the the number of values read
     */
    template<class T> void sethist(T* h, const std::vector<double> &values, size_t &index){
        const int n = h->GetSize();
        //h->Set(n, &values[index]);
        //index += (size_t)n;
        for(int i=0; i<n; ++i){
            h->SetBinContent(i,values[index++]);
        }
    }

    /**
     * Add the values of a histogram to the std::vector
     */
    template<class T> void addhist(T* h, std::vector<double> &values){
        const int n = h->GetSize();
        const double *array = h->GetArray();
        values.insert(values.end(), array, array + n);
    }
}

DeltaTHists::DeltaTHists(int nbins, double mindt, double maxdt, const std::string &name):deleteHists(false){
    boost::format dt_hist_name("dT_svd%1%_rbin%2%_fit_q%3%_e%4%_%5%");
    boost::format dt_hist_title("#Deltat, SVD%1%, %5%, rbin %2%, q=%3%1, #eta=%4%1");
    const std::string pmn("pm");
    const std::string pmt("+-");

    if(name.empty()){
        //Do not save the histograms
        TH1::AddDirectory(false);
        deleteHists = true;
    }

    for(int svd=0; svd<NSVD; ++svd){
        for(int rbin=0; rbin<NRBIN; ++rbin){
            for(int q=0; q<2; ++q){
                for(int eta=0; eta<2; ++eta){
                    const std::string hname = (dt_hist_name % (svd+1) % rbin % pmn[q] % pmn[eta] % name).str();
                    const std::string htitle = (dt_hist_title % (svd+1) % rbin % pmt[q] % pmt[eta] % name).str();
                    hists[index(svd, rbin, q, eta)] = new TH1D(hname.c_str(),htitle.c_str(), nbins, mindt, maxdt);
                }
            }
        }
    }

    yields = new TH2D(("nevents_" + name).c_str(), "#Events;SVD;rbin", NSVD, 0, NSVD, NRBIN, 0, NRBIN);

    TH1::AddDirectory(true);
}

DeltaTHists::~DeltaTHists(){
    if(deleteHists){
        delete yields;
        for(int i=0; i<NHISTS; ++i){
            delete hists[i];
        }
    }
}

void DeltaTHists::finalize(std::function<double (int, int)> pdf_yield){
    const double binsize = hists[0]->GetXaxis()->GetBinWidth(1);
    for(int svd=0; svd<NSVD; ++svd){
        for(int rbin=0; rbin<NRBIN; ++rbin){
            const double nevents = yield(svd,rbin);
            if(nevents==0) continue;
            for(int q=0; q<2; ++q){
                for(int eta=0; eta<2; ++eta){
                    TH1D* h = hists[index(svd,rbin,q,eta)];
                    h->Scale(pdf_yield(svd,rbin) * binsize / nevents);
                }
            }
        }
    }
}

void DeltaTHists::recieve(const std::vector<double> &values){
    size_t index(0);
    sethist(yields, values, index);
    for(int i=0; i<NHISTS; ++i){
        sethist(hists[i], values, index);
    }
}

std::vector<double> DeltaTHists::send() {
    std::vector<double> values;
    values.reserve(yields->GetSize() + NHISTS*hists[0]->GetSize());
    addhist(yields, values);
    for(int i=0; i<NHISTS; ++i){
        addhist(hists[i], values);
    }
    return values;
}
