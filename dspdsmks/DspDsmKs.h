#ifndef MPIFitter_DspDsmKsPDF_h
#define MPIFitter_DspDsmKsPDF_h

#include <functional>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
//#include "func.h"

#include <TChain.h>
#include <Functions.h>
#include "Range.h"
#include "Event.h"
#include "Signal.h"
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
        CMP_NONE        = 0,
        CMP_signal      = 1<<0,
        CMP_mixed       = 1<<1,
        CMP_charged     = 1<<2,
        CMP_dt_signal   = 1<<3,
        CMP_dt_mixed    = 1<<4,
        CMP_dt_charged  = 1<<5,
        CMP_dt_all      = CMP_dt_signal | CMP_dt_mixed | CMP_dt_charged,
        CMP_nodt_all    = CMP_signal | CMP_mixed | CMP_charged,
        CMP_all         = CMP_dt_all | CMP_nodt_all
    };

    static EnabledComponents getComponents(const std::vector<std::string> &components){
        EnabledComponents result = CMP_NONE;
#define DspDsmKsPDF__checkComponent(name) if(component == #name) result = (EnabledComponents) (result | CMP_##name)
        BOOST_FOREACH(std::string component, components){
            boost::to_lower(component);
            boost::trim(component);
            DspDsmKsPDF__checkComponent(signal);
            DspDsmKsPDF__checkComponent(mixed);
            DspDsmKsPDF__checkComponent(charged);
            DspDsmKsPDF__checkComponent(dt_signal);
            DspDsmKsPDF__checkComponent(dt_mixed);
            DspDsmKsPDF__checkComponent(dt_charged);
            DspDsmKsPDF__checkComponent(dt_all);
            DspDsmKsPDF__checkComponent(nodt_all);
            DspDsmKsPDF__checkComponent(all);
        }
#undef DspDsmKsPDF__checkComponent
        return result;
    }

    DspDsmKsPDF(Range range_mBC, Range range_dE, Range range_dT,
            const std::vector<std::string> filenames, const std::string &bestB, EnabledComponents components=CMP_all, int maxPrintOrder = 0):
        svdVs(Component::BOTH), maxPrintOrder(maxPrintOrder), nCalls(0), filenames(filenames), bestBSelection(bestB),
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
            components.push_back(new SignalPDF(range_mBC, range_dE, range_dT, cmp & CMP_dt_signal));
        }
        if(cmp & CMP_mixed){
            components.push_back(new MixedPDF(range_mBC, range_dE, range_dT, cmp & CMP_dt_mixed));
        }
        if(cmp & CMP_charged){
            components.push_back(new ChargedPDF(range_mBC, range_dE, range_dT, cmp & CMP_dt_charged));
        }
    }

    void setSVD(Component::EnabledSVD svd){
        svdVs = svd;
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

        for( unsigned int i=0; i<data.size(); i++ ){
            const double pdf = PDF(data[i], par);
            if( 0.0 < pdf ) log_pdf += log( pdf ) ;
        }
        return log_pdf;
    }

    /** Return the yield of the pdf given the set of parameters */
    double yield(const std::vector<double> &par) const {
        double yield(0);
        BOOST_FOREACH(Component* component, components){
            yield += component->get_yield(par,svdVs);
        }
        return yield;
    }

    /** finalize the event after all processes are collected */
    double finalize(const std::vector<double> &par, double value) const {
        const double logL = -2.0*(value-yield(par));

        //Determine wether to show value
        nCalls++;
        const int order = (nCalls == 0) ? 1 : std::max(std::min((int)std::log10(nCalls), maxPrintOrder), 0);
        const int interval = static_cast<int>(std::pow(10., order));
        if(nCalls % interval == 0){
            std::cout << "call #" << std::setw(5) << std::left << nCalls << ": ";
            std::cout << "-2logL =" << std::setw(18) << std::setprecision(10) <<
                std::scientific << std::right << logL << std::endl;
        }
        return logL;
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
        const unsigned int entries = chain->GetEntries();
        const unsigned int interval = std::ceil(1.0 * entries / size);
        const unsigned int start = process*interval;
        const unsigned int end = std::min(entries,start+interval);

        data.clear();
        data.reserve(end-start+1);
        Event event;
        event.setBranches(chain,bestBSelection);
        for(unsigned int i=start; i<end; ++i){
            if(!chain->GetEntry(i)) {
                std::cerr << "ARRRRRRRRR: There be a problem readin in event " << i << std::endl;
                std::exit(5);
            }
            if(!range_mBC(event.Mbc) || !range_dE(event.dE)) continue;
            event.calculateValues();
            data.push_back(event);
        }

        std::cout << "Aye, process " << process << " fully loaded " << data.size()
            << " events and is ready for pillaging" << std::endl;
        delete chain;
    }

    const std::vector<Event>& getData() const { return data; }

    protected:

    /** Which components to include into the pdf */
    //EnabledComponents components;
    Component::EnabledSVD svdVs;
    /** Max order to print out log2L */
    int maxPrintOrder;
    /** Number of calls */
    mutable int nCalls;
    /** vector containing the data values */
    std::vector<Event> data;
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
