#ifndef MPIFitter_DspDsmKsPDF_h
#define MPIFitter_DspDsmKsPDF_h

#include <functional>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <boost/random/random_device.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/poisson_distribution.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <algorithm>

#include <TChain.h>
#include <Functions.h>
#include <TH1D.h>
#include "Range.h"
#include "Event.h"
#include "Signal.h"
#include "Misrecon.h"
#include "Mixed.h"
#include "Charged.h"
#include "progress.h"

/** DspDsmKs PDF function.
 *
 * This is a simple 1D, 2components fit from jeremys tutorial
 */
struct DspDsmKsPDF {

    /** Define how different pieces of the pdf should be reduced to one value
     *
     * useful possibilities include:
     * std::plus<double> -> add the value from all processes
     * std::multiplies<double> -> multiply the value from all processes
     */
    typedef std::plus<double> operator_type;

    /** Define the increase of the objective function value which is eqvivalent to 1 sigma */
    constexpr static double error_def = 1.0;

    /** Define all possible components to be able to enable/disable them individually */
    enum EnabledComponents {
        CMP_NONE     = 0,
        CMP_signal   = 1<<0,
        CMP_misrecon = 1<<1,
        CMP_mixed    = 1<<2,
        CMP_charged  = 1<<3,
        CMP_deltat   = 1<<4,
        CMP_all      = CMP_signal | CMP_misrecon | CMP_mixed | CMP_charged | CMP_deltat
    };

    enum PlotFlags {
        PLT_NONE   = 0,
        PLT_SVD1   = 1<<0,
        PLT_SVD2   = 1<<1,
        PLT_DT_Q   = 1<<2,
        PLT_DT_E   = 1<<3,
        PLT_DT_QE  = 1<<4,
        PLT_DT     = PLT_DT_Q | PLT_DT_E | PLT_DT_QE,
        PLT_MAX    = 1<<5,
        PLT_MAXDT  = 1<<6
    };

    static EnabledComponents getComponents(const std::string &components){
        EnabledComponents result = CMP_NONE;
#define DspDsmKsPDF__checkComponent(name) if(component == #name) result = (EnabledComponents) (result | CMP_##name)
        boost::char_separator<char> sep(", ");
        boost::tokenizer<boost::char_separator<char> > tokens(components,sep);
        for(std::string component: tokens){
            boost::to_lower(component);
            boost::trim(component);
            DspDsmKsPDF__checkComponent(signal);
            else DspDsmKsPDF__checkComponent(misrecon);
            else DspDsmKsPDF__checkComponent(mixed);
            else DspDsmKsPDF__checkComponent(charged);
            else DspDsmKsPDF__checkComponent(deltat);
            else DspDsmKsPDF__checkComponent(all);
            else throw std::invalid_argument("Unknown component: '" + component + "'");
        }
#undef DspDsmKsPDF__checkComponent
        return result;
    }

    DspDsmKsPDF(int maxPrintOrder = 2):
        maxPrintOrder(maxPrintOrder), nCalls(0), minLogL(std::numeric_limits<double>::infinity()), maxLogL(-std::numeric_limits<double>::infinity()), bestBSelection("bestLHsig"),
        range_mBC(5.24,5.30), range_dE(-0.15,0.1), range_dT(Belle::dt_resol_global::dt_llmt, Belle::dt_resol_global::dt_ulmt)
    {}

    ~DspDsmKsPDF(){
        for(Component* component: components){
            delete component;
        }
    }

    void setComponents(EnabledComponents cmp=CMP_all){
        enabledComponents = cmp;
        for(Component* component: components){
            delete component;
        }
        components.clear();
        if(cmp & CMP_signal){
            components.push_back(new SignalPDF(range_mBC, range_dE, range_dT, cmp & CMP_deltat));
        }
        if(cmp & CMP_misrecon){
            components.push_back(new MisreconPDF(range_mBC, range_dE, range_dT, cmp & CMP_deltat));
        }
        if(cmp & CMP_mixed){
            components.push_back(new MixedPDF(range_mBC, range_dE, range_dT, cmp & CMP_deltat));
        }
        if(cmp & CMP_charged){
            components.push_back(new ChargedPDF(range_mBC, range_dE, range_dT, cmp & CMP_deltat));
        }
    }

    void setOptions(int options) {
        setComponents((EnabledComponents) options);
    }

    /** Return the pdf value for a given paramater set and event */
    double PDF(const Event& e, const std::vector<double> &par) const {
        long double result(0);
        for(Component* component: components){
            result += (*component)(e, par);
        }
        return result;
    }

    /** Return the pdf for all events */
    double operator()(const std::vector<double> &par) const {
        double log_pdf(0.0);
        for(int svd=0; svd<2; ++svd){
            for( unsigned int i=0; i<data[svd].size(); i++ ){
                const double pdf = PDF(data[svd][i], par);
                if( 0.0 < pdf ) log_pdf += log( pdf ) ;
            }
        }
        return log_pdf;
    }

    /** Return the yield of the pdf given the set of parameters */
    double get_yield(const std::vector<double> &par, int svdVs=Component::BOTH, int rbin=-1) const {
        double yield(0);
        for(Component* component: components){
            yield += component->get_yield(par,(Component::EnabledSVD)svdVs, rbin);
        }
        return yield;
    }

    double get_deltaT(const Event &e, const std::vector<double> &par) const {
        long double deltaT(0);
        for(Component* component: components){
            deltaT += component->get_deltaT(e,par, true);
        }
        return deltaT;
    }

    std::vector<double> get_rbinFractions(const std::vector<double> &par, int svd){
        std::vector<double> result(7);
        std::vector<double> fractions(7);
        double yield(0);
        for(Component* component: components){
            const double y = component->get_yield(par, (Component::EnabledSVD)svd);
            component->get_rbinFractions(par, fractions, (Component::EnabledSVD)svd, y);
            for(int i=0; i<7; ++i){ result[i] += fractions[i]; }
            yield += y;
        }
        for(int i=0; i<7; ++i){ result[i] /= yield; }
        return result;
    }


    double get_cosTheta(const Event &e, const std::vector<double> &par, int svdVs) const {
        long double yield(0);
        long double sum(0);
        for(Component* component: components){
            const double y = component->get_yield(par,(Component::EnabledSVD)svdVs);
            const double c = component->get_cosTheta(e);
            sum += y*c;
            yield += y;
        }
        return sum/yield;
    }

    /** finalize the event after all processes are collected */
    double finalize(const std::vector<double> &par, double value) const {
        static boost::format output("call #%-5d: -2logL =%18.10g     (min =%18.10g, max =%18.10g)\n");
        const double logL = -2.0*(value-get_yield(par));
        if(logL<minLogL) minLogL = logL;
        if(logL>maxLogL) maxLogL = logL;

        //Determine wether to show value
        nCalls++;
        const int order = (nCalls == 0) ? 1 : std::max(std::min((int)std::log10(nCalls), maxPrintOrder), 0);
        const int interval = static_cast<int>(std::pow(10., order));
        if(nCalls % interval == 0){
            std::cout << output % nCalls % logL % minLogL % maxLogL;
        }
        return logL;
    }

    std::vector<double> plotM(int flag, const std::vector<double> values, const std::vector<double> &par, int rank){
        int svdVs=0;
        std::vector<double> results;
        if(flag & PLT_SVD1) svdVs |= Component::SVD1;
        if(flag & PLT_SVD2) svdVs |= Component::SVD2;

        if(flag & PLT_DT){
            TH1D *h_dt[7][2][2];
            for(int rbin=0; rbin<7; ++rbin){
                for(int i=0; i<2; ++i){
                    for(int j=0; j<2; ++j){
                        const std::string name = (boost::format("dt_%d_%d_%d") % rbin % i % j).str();
                        h_dt[rbin][i][j] = new TH1D(name.c_str(),"",(int)values[0],values[1],values[2]);
                    }
                }
            }
            for(int svd=0; svd<2; ++svd){
                if((svd==0) && !(flag & PLT_SVD1)) continue;
                if((svd==1) && !(flag & PLT_SVD2)) continue;
                ProgressBar pbar(h_dt[0][0][0]->GetNbinsX()*data[svd].size());
                if(rank==0 && data[svd].size()==0) std::cout << ++pbar;
                for(Event e: data[svd]){
                    for(int ix=0; ix<h_dt[0][0][0]->GetNbinsX(); ++ix){
                        if(rank==0) std::cout << ++pbar;
                        e.deltaT = h_dt[0][0][0]->GetBinCenter(ix+1);
                        e.reset();
                        for(int q=0; q<2; ++q){
                            e.tag_q = 2*q-1;
                            for(int eta=0; eta<2; ++eta){
                                e.eta = 2*eta-1;
                                double pdf_value = get_deltaT(e, par);
                                h_dt[e.rbin][0][q]->Fill(e.deltaT, pdf_value);
                                h_dt[e.rbin][1][(q!=eta)?0:1]->Fill(e.deltaT, pdf_value);
                            }
                        }
                    }
                }
            }
            for(int rbin=0; rbin<7; ++rbin){
                for(int i=0; i<2; ++i){
                    double h_scale = h_dt[rbin][i][0]->GetXaxis()->GetBinWidth(1);
                    for(int j=0; j<2; ++j){
                        TH1D* h = h_dt[rbin][i][j];
                        h->Scale(h_scale);
                        const double* array = h->GetArray();
                        //Stupid underflow bin is array[0], so we need to start at element 1
                        results.insert(results.end(), array + 1, array + h->GetNbinsX() + 1);
                        delete h;
                    }
                }
            }
        }
        return results;
    }

    double plot(int flag, const std::vector<double> values, const std::vector<double> &par){
        int svdVs=0;
        if(flag & PLT_SVD1) svdVs |= Component::SVD1;
        if(flag & PLT_SVD2) svdVs |= Component::SVD2;

        if(flag & PLT_DT){
            long double pdf(0.0);
            int nEvents(0);
            for(int svd=0; svd<2; ++svd){
                if((svd==0) && !(flag & PLT_SVD1)) continue;
                if((svd==1) && !(flag & PLT_SVD2)) continue;
                for(unsigned int i=0; i<data[svd].size(); i++ ){
                    Event e = data[svd][i];
                    e.deltaT = values[0];
                    for(int i=-1; i<=2; i+=2){
                        if(flag & PLT_DT_Q){
                            e.tag_q = (int) values[1];
                            e.eta = i;
                        }else if(flag & PLT_DT_E){
                            e.tag_q = i;
                            e.eta = (int) values[1];
                        }else if(flag & PLT_DT_QE){
                            e.tag_q = i;
                            e.eta = (int) (i*values[1]);
                        }
                        pdf += get_deltaT(e, par);
                    }
                    ++nEvents;
                }
            }
            const double yield = get_yield(par,svdVs);
            return pdf*yield/nEvents;
        }
        if(flag & (PLT_MAX | PLT_MAXDT)){
            double pdf(0.0);
            int oldComponents = enabledComponents;
            int svd = (flag & PLT_SVD1)?0:1;
            Event e;
            e.svdVs = svd;
            e.benergy = values[0];
            if(flag & PLT_MAX){
                if(oldComponents & CMP_deltat) setOptions(oldComponents ^ CMP_deltat);
                const double &mbc_start = values[1];
                const double &mbc_end   = values[2];
                const double &mbc_step  = values[3];
                const double &de_start  = values[4];
                const double &de_end    = values[5];
                const double &de_step   = values[6];
                for(e.Mbc=mbc_start; e.Mbc <= mbc_end; e.Mbc+=mbc_step){
                    for(e.dE=de_start; e.dE <= de_end; e.dE+=de_step){
                        pdf = std::max(pdf, PDF(e,par));
                    }
                }
            }else{
                pdf = 1.0;
            }
            if(oldComponents & CMP_deltat){
                const double &dt_start  = values[7];
                const double &dt_end    = values[8];
                const double &dt_step   = values[9];
                double dtpdf(0);
                std::string name = "dt_max_svd";
                name += svd?"2":"1";
                for(double dt=dt_start; dt <= dt_end; dt+=dt_step){
                    for(unsigned int i=0; i<std::min((size_t)100u,data[svd].size()); i++){
                        e = data[svd][i];
                        e.deltaT = dt;
                        for(e.tag_q=-1; e.tag_q<=2; e.tag_q+=2){
                            for(e.eta=-1; e.eta<=2; e.eta+=2){
                                dtpdf = std::max(dtpdf,get_deltaT(e,par));
                            }
                        }
                    }
                }
                pdf *= dtpdf;
                if(flag & PLT_MAX) setOptions(oldComponents);
            }
            return pdf;
        }
        return 0;
    }

    void selectToyMC(TTree* output, const std::vector<double> &par, int seed=0){
        boost::random::mt19937 random_generator(seed);
        if(seed == 0){
            boost::random::random_device rseed;
            random_generator.seed(rseed);
        }
        Event e;
        e.createBranches(output, bestBSelection);
        typedef boost::random::uniform_int_distribution<> uniform_int;
        typedef boost::variate_generator<boost::random::mt19937&, uniform_int> int_variate;

        for(int svd=0; svd<2; ++svd){
            int_variate  random_event(random_generator, uniform_int(0, data[svd].size()-1));
            double yield = get_yield(par, svd==0?Component::SVD1:Component::SVD2);
            std::cout << "yield for SVD " << svd << ": " << yield << std::endl;
            if(yield<=0) continue;
            boost::random::poisson_distribution<> poisson(yield);
            int nEvents = poisson(random_generator);
            if(nEvents<=0) continue;
            std::cout << "Selecting " << nEvents << " Events for SVD" << (svd+1) << std::endl;
            for(int i=0; i<nEvents; ++i){
                e = data[svd][random_event()];
                output->Fill();
            }
        }
    }

    void generateToyMC(TTree* output, const std::vector<double> &par, double maxval[2], std::vector<std::string> &templates, int seed=0, bool fullgsim=false){
        boost::random::mt19937 random_generator(seed);
        if(seed == 0){
            boost::random::random_device rseed;
            random_generator.seed(rseed);
        }

        bool gsim = false;
        std::vector<Event> gsim_data[2];
        if(!templates.empty()){
            std::cout << "Generating from gsim" << std::endl;
            loadEvents(gsim_data, templates, 0, 1, false);
            gsim = true;
        }

        Event e;
        e.createBranches(output, bestBSelection);
        e.flag = 0;
        typedef boost::random::uniform_int_distribution<> uniform_int;
        typedef boost::random::uniform_real_distribution<> uniform_real;
        typedef boost::variate_generator<boost::random::mt19937&, uniform_int> int_variate;
        typedef boost::variate_generator<boost::random::mt19937&, uniform_real > real_variate;
        typedef boost::variate_generator<boost::random::mt19937&, boost::random::discrete_distribution<> > rbin_variate;
        for(int svd=0; svd<2; ++svd){
            if(gsim && gsim_data[svd].empty()){
                std::cerr << "No template data for SVD" << (svd+1) << std::endl;
                std::exit(5);
            }
            int_variate  random_event(random_generator, uniform_int(0, data[svd].size()-1));
            int_variate  random_gsim(random_generator, uniform_int(0, std::max(0,(int)(gsim_data[svd].size())-1)));
            int_variate  random_flavour(random_generator, uniform_int(0, 1));
            real_variate random_pdf(random_generator, uniform_real(0, maxval[svd]));
            real_variate random_mBC(random_generator, uniform_real(range_mBC.vmin, range_mBC.vmax));
            real_variate random_dE(random_generator, uniform_real(range_dE.vmin, range_dE.vmax));
            real_variate random_dT(random_generator, uniform_real(range_dT.vmin, range_dT.vmax));
            real_variate random_cos(random_generator, uniform_real(-1,1));
            real_variate random_01(random_generator, uniform_real(0,1));

            int nEvents = (int)round(get_yield(par, svd==0?Component::SVD1:Component::SVD2));
            std::cout << "yield for SVD " << svd << ": " << nEvents << std::endl;
            if(nEvents<=0) continue;
            boost::random::poisson_distribution<> poisson(nEvents);
            nEvents = poisson(random_generator);
            if(nEvents==0) continue;

            std::vector<Event> rbinpool[7];
            for(const Event& e: data[svd]){
                rbinpool[e.rbin].push_back(e);
            }
            int_variate random_rbinevent[7] = {
                int_variate(random_generator, uniform_int(0, rbinpool[0].size()-1)),
                int_variate(random_generator, uniform_int(0, rbinpool[1].size()-1)),
                int_variate(random_generator, uniform_int(0, rbinpool[2].size()-1)),
                int_variate(random_generator, uniform_int(0, rbinpool[3].size()-1)),
                int_variate(random_generator, uniform_int(0, rbinpool[4].size()-1)),
                int_variate(random_generator, uniform_int(0, rbinpool[5].size()-1)),
                int_variate(random_generator, uniform_int(0, rbinpool[6].size()-1)),
            };
            std::vector<Event> rbingsim[7];
            for(const Event& e: gsim_data[svd]){
                rbingsim[e.rbin].push_back(e);
            }
            int_variate random_gsimevent[7] = {
                int_variate(random_generator, uniform_int(0, rbingsim[0].size()-1)),
                int_variate(random_generator, uniform_int(0, rbingsim[1].size()-1)),
                int_variate(random_generator, uniform_int(0, rbingsim[2].size()-1)),
                int_variate(random_generator, uniform_int(0, rbingsim[3].size()-1)),
                int_variate(random_generator, uniform_int(0, rbingsim[4].size()-1)),
                int_variate(random_generator, uniform_int(0, rbingsim[5].size()-1)),
                int_variate(random_generator, uniform_int(0, rbingsim[6].size()-1)),
            };

            std::vector<double> rbinFractions = get_rbinFractions(par, svd==0?Component::SVD1:Component::SVD2);
            for(int i=0; i<7; ++i){
                std::cout << "SVD" << (svd+1) << ", rbin " << i
                    << ": fraction=" << rbinFractions[i] << ", events=" << rbinpool[i].size() << std::endl;
            }
            rbin_variate random_rbin(random_generator, boost::random::discrete_distribution<>(rbinFractions.begin(), rbinFractions.end()));

            e.svdVs = svd;
            e.isMC = data[svd][0].isMC;
            ProgressBar pbar(nEvents);
            std::cout << "Generating " << nEvents << " Events for SVD" << (svd+1) << ":" << pbar;
            double min_distance(std::numeric_limits<double>::infinity());
            for(int i=0; i<nEvents; ++i){
                //Take vertex stuff from data
                do {
                    e.rbin = random_rbin();
                    std::vector<Event> &pool = rbinpool[e.rbin];
                    int_variate &pool_event = random_rbinevent[e.rbin];
                    //e.cosTheta = data[svd][random_event()].cosTheta;
                    //e.tag_r    = pool[pool_event()].tag_r;
                    e.tag_zerr = pool[pool_event()].tag_zerr;
                    e.tag_chi2 = pool[pool_event()].tag_chi2;
                    e.tag_ndf  = pool[pool_event()].tag_ndf;
                    e.tag_isL  = pool[pool_event()].tag_isL;
                    e.vtx_zerr = pool[pool_event()].vtx_zerr;
                    e.vtx_chi2 = pool[pool_event()].vtx_chi2;
                    e.vtx_ndf  = pool[pool_event()].vtx_ndf;
                    e.tag_ntrk = (e.tag_ndf+2)/2;
                    e.vtx_ntrk = (e.vtx_ndf+2)/2;
                }while(!e.calculateValues(true, true));

                //Take Mbc/dE stuff from template if available
                if(gsim){
                    Event &e2 = rbingsim[e.rbin][random_gsimevent[e.rbin]()];
                    //Event &e2 = gsim_data[svd][random_gsim()];
                    e.benergy = e2.benergy;
                    e.expNo = e2.expNo;
                    e.Mbc = e2.Mbc;
                    e.dE = e2.dE;
                    //When not generating signal take the dt parameters from data too.
                    if(fullgsim && !(enabledComponents & (CMP_signal | CMP_misrecon))){
                        e.eta = e2.eta;
                        e.tag_q = e2.tag_q;
                        e.deltaT = e2.deltaT;
                        assert(e.calculateValues(true, true));
                        //Everything set, go to next event
                        std::cout << ++pbar;
                        output->Fill();
                        continue;
                    }
                }

                while(true){
                    double calc_pdf(-1);
                    if(!gsim){
                        Event &e2 = data[svd][random_event()];
                        e.benergy = e2.benergy;
                        e.expNo = e2.expNo;
                        //Make sure Mbc is in a valid range
                        do{
                            e.Mbc = random_mBC();
                            e.dE = random_dE();
                        }while(e.Mbc>e.benergy);
                    }

                    //No we get a cosTheta
                    do { e.cosTheta = random_cos(); }
                    while(get_cosTheta(e,par,svd==0?Component::SVD1:Component::SVD2)<random_01());

                    //And now we make deltaT
                    if(enabledComponents & CMP_deltat){
                        e.deltaT   = random_dT();
                        e.tag_q    = (random_flavour()*2)-1;
                        e.eta      = (random_flavour()*2)-1;
                        assert(e.calculateValues(true, true));
                        e.reset();
                    }
                    //If gsim and there is no deltaT, nothing else to do, just grab a random event
                    if(gsim && !(enabledComponents & CMP_deltat)) break;
                    //otherwise get pdf/dt value and check against random value
                    calc_pdf = gsim?get_deltaT(e,par):PDF(e,par);

                    if(calc_pdf>maxval[svd]){
                        throw std::runtime_error(
                                (boost::format("PDF larger than maximimum value: %.10e > %.10e") % calc_pdf % maxval[svd]).str());
                    }
                    min_distance = std::min((maxval[svd]-calc_pdf)/maxval[svd], min_distance);

                    if(random_pdf()<calc_pdf) break;
                }

                std::cout << ++pbar;
                output->Fill();
            }
            std::cout << "Min distance to max value: " << min_distance << std::endl;
        }
    }

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
    void load(int process=0, int size=1, bool qualityCuts=true) {
        loadEvents(data, filenames, process, size, qualityCuts);
        std::sort(data[0].begin(),data[0].end());
        std::sort(data[1].begin(),data[1].end());
        std::cout << "Aye, process " << process << " fully loaded (" << data[0].size() << ", " << data[1].size()
            << ") events and is ready for pillaging" << std::endl;
    }

    const std::vector<Event>& getData(int svd) const { return data[svd]; }
    size_t size(int svd) const { return data[svd].size(); }

    const Range getRange_mBC() const { return range_mBC; }
    const Range getRange_dE() const { return range_dE; }
    const Range getRange_dT() const { return range_dT; }

    Range& getRange_mBC() { return range_mBC; }
    Range& getRange_dE() { return range_dE; }
    Range& getRange_dT() { return range_dT; }

    std::string& getBestB() { return bestBSelection; }
    std::vector<std::string>& getFiles(){ return filenames; }

    protected:

    void loadEvents(std::vector<Event> *data, const std::vector<std::string>& filenames, int process, int size, bool qualityCuts){
        TChain* chain = new TChain("B0");
        for(const std::string& filename: filenames){
            std::cout << "Using " << filename << " to read events" << std::endl;
            chain->AddFile(filename.c_str(),-1);
        }
        data[0].clear();
        data[1].clear();
        Event event;
        event.setBranches(chain, bestBSelection);
        const unsigned int entries = chain->GetEntries();
        for(unsigned int i=process; i<entries; i+=size){
            if(!chain->GetEntry(i)) {
                std::cerr << "ARRRRRRRRR: There be a problem readin in event " << i << std::endl;
                std::exit(5);
            }
            if(!event.calculateValues(qualityCuts)) continue;
            if(!range_mBC(event.Mbc) || !range_dE(event.dE) || !range_dT(event.deltaT)) continue;
            data[event.svdVs].push_back(event);
        }
        delete chain;
    }


    /** Max order to print out log2L */
    int maxPrintOrder;
    /** Number of calls */
    mutable int nCalls;
    mutable double minLogL;
    mutable double maxLogL;
    /** vector containing the data values */
    std::vector<Event> data[2];
    /** filename from which the data should be read */
    std::vector<std::string> filenames;
    /** name of the bestBSelection to be used */
    std::string bestBSelection;

    /** Range of the mBC fit */
    Range range_mBC;
    /** Range of the dE fit */
    Range range_dE;
    /** Range of the dT fit */
    Range range_dT;

    /** PDF function components */
    mutable std::vector<Component*> components;
    EnabledComponents enabledComponents;
};

#endif //MPIFitter_DspDsmKsPDF_h
