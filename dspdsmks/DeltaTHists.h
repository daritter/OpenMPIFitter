#include <vector>
#include <string>
#include <functional>

#include <TH1D.h>
#include <TH2D.h>

#include "Event.h"

class DeltaTHists {
    public:
        enum {
            NSVD = 2,
            NRBIN = 7,
            NHISTS_PER_RBIN = 4,
            NHISTS_PER_SVD = NRBIN*NHISTS_PER_RBIN,
            NHISTS = NSVD*NHISTS_PER_SVD
        };

        DeltaTHists(int nbins, double mindt, double maxdt, const std::string &name="");

        ~DeltaTHists();

        void finalize();

        void recieve(const std::vector<double> &values);

        std::vector<double> send();

        void fill(Event e, std::function<double (const Event&)> pdf){
            yields->Fill(e.svdVs, e.rbin);
            for(int ix=0; ix<hists[0]->GetNbinsX(); ++ix){
                e.deltaT = hists[0]->GetBinCenter(ix+1);
                for(int q=0; q<2; ++q){
                    e.tag_q = 2*q-1;
                    for(int eta=0; eta<2; ++eta){
                        e.eta = 2*eta-1;
                        hists[index(e.svdVs, e.rbin, q, eta)]->Fill(e.deltaT, pdf(e));
                    }
                }
            }
        }

    protected:
        int index(int svd, int rbin, int q, int eta){
            return svd*NHISTS_PER_SVD + rbin*NHISTS_PER_RBIN + q*2 + eta;
        }

        bool deleteHists;// = false;
        TH2D* yields;// = 0;
        TH1D* hists[NHISTS];// = {0};
};
