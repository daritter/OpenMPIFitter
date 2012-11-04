#ifndef MPIFitter_DspDsmKsPDF_h
#define MPIFitter_DspDsmKsPDF_h

#include <functional>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <algorithm>
//#include "func.h"

#include <TChain.h>
#include <Functions.h>
#include "Range.h"
#include "Event.h"
#include "Signal.h"
#include "Misrecon.h"
#include "Mixed.h"
#include "Charged.h"

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
    static const double error_def = 1.0;

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
        PLT_DT_QE  = 1<<4
    };

    static EnabledComponents getComponents(const std::string &components){
        EnabledComponents result = CMP_NONE;
#define DspDsmKsPDF__checkComponent(name) if(component == #name) result = (EnabledComponents) (result | CMP_##name)
        boost::char_separator<char> sep(", ");
        boost::tokenizer<boost::char_separator<char> > tokens(components,sep);
        BOOST_FOREACH(std::string component, tokens){
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

    DspDsmKsPDF(Range range_mBC, Range range_dE, Range range_dT,
            const std::vector<std::string> filenames, const std::string &bestB, EnabledComponents components=CMP_all, int maxPrintOrder = 0):
        maxPrintOrder(maxPrintOrder), nCalls(0), filenames(filenames), bestBSelection(bestB),
        range_mBC(range_mBC), range_dE(range_dE), range_dT(range_dT)
    {
        setComponents(components);
    }

    ~DspDsmKsPDF(){
        BOOST_FOREACH(Component* component, components){
            delete component;
        }
    }

    void setComponents(EnabledComponents cmp=CMP_all){
        BOOST_FOREACH(Component* component, components){
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
        BOOST_FOREACH(Component* component, components){
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
    double get_yield(const std::vector<double> &par, int svdVs=Component::BOTH) const {
        double yield(0);
        BOOST_FOREACH(Component* component, components){
            yield += component->get_yield(par,(Component::EnabledSVD)svdVs);
        }
        return yield;
    }

    double get_deltaT(const Event &e, const std::vector<double> &par) const {
        long double deltaT(0);
        BOOST_FOREACH(Component* component, components){
            deltaT += component->get_deltaT(e,par, true);
        }
        return deltaT;
    }

    /** finalize the event after all processes are collected */
    double finalize(const std::vector<double> &par, double value) const {
        static boost::format output("call #%-5d: -2logL =%18.10g\n");
        const double logL = -2.0*(value-get_yield(par));

        //Determine wether to show value
        nCalls++;
        const int order = (nCalls == 0) ? 1 : std::max(std::min((int)std::log10(nCalls), maxPrintOrder), 0);
        const int interval = static_cast<int>(std::pow(10., order));
        if(nCalls % interval == 0){
            std::cout << output % nCalls % logL;
        }
        return logL;
    }

    double plot(int flag, const std::vector<double> values, const std::vector<double> &par){
        long double pdf(0.0);
        int svdVs=0;
        if(flag & PLT_SVD1) svdVs |= Component::SVD1;
        if(flag & PLT_SVD2) svdVs |= Component::SVD2;
        int nEvents(0);
        for(int svd=0; svd<2; ++svd){
            if((svd==0) && !(flag & PLT_SVD1)) continue;
            if((svd==1) && !(flag & PLT_SVD2)) continue;
            for(unsigned int i=0; i<data[svd].size(); i++ ){
                Event e = data[svd][i];
                e.deltaT = values[0];
                for(int i=-1; i<=2; i+=2){
                    if(flag & PLT_DT_Q){
                        e.tag_q = values[1];
                        e.eta = i;
                    }else if(flag & PLT_DT_E){
                        e.tag_q = i;
                        e.eta = values[1];
                    }else if(flag & PLT_DT_QE){
                        e.tag_q = i;
                        e.eta = i*values[1];
                    }
                    pdf += get_deltaT(e, par);
                }
                ++nEvents;
            }
        }
        const double yield = get_yield(par,svdVs);
        return pdf*yield/nEvents;
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
    void load(int process, int size) {
        TChain* chain = new TChain("B0");

        BOOST_FOREACH(const std::string& filename, filenames){
            chain->AddFile(filename.c_str(),-1);
        }
        bool readStriped = true;
        const unsigned int entries = chain->GetEntries();
        unsigned int start=process;
        unsigned int end=entries;
        unsigned int stride=size;
        if(!readStriped){
            const unsigned int interval = std::ceil(1.0 * entries / size);
            start = process*interval;
            end = std::min(entries,start+interval);
            stride = 1;
        }

        data[0].clear();
        data[1].clear();
        //data.reserve(end-start+1);
        Event event;
        event.setBranches(chain,bestBSelection);
        for(unsigned int i=start; i<end; i+=stride){
            if(!chain->GetEntry(i)) {
                std::cerr << "ARRRRRRRRR: There be a problem readin in event " << i << std::endl;
                std::exit(5);
            }
            event.calculateValues();
            if(!range_mBC(event.Mbc) || !range_dE(event.dE) || !range_dT(event.deltaT)) continue;
            if(event.flag!=0 || event.tag_q==0) continue;
            data[event.svdVs].push_back(event);
        }

        std::sort(data[0].begin(),data[0].end());
        std::sort(data[1].begin(),data[1].end());

        std::cout << "Aye, process " << process << " fully loaded (" << data[0].size() << ", " << data[1].size()
            << ") events and is ready for pillaging" << std::endl;
        delete chain;
    }

    const std::vector<Event>& getData(int svd) const { return data[svd]; }


    const Range getRange_mBC() const { return range_mBC; }
    const Range getRange_dE() const { return range_dE; }
    const Range getRange_dT() const { return range_dT; }

    protected:

    /** Max order to print out log2L */
    int maxPrintOrder;
    /** Number of calls */
    mutable int nCalls;
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
};

#endif //MPIFitter_DspDsmKsPDF_h
