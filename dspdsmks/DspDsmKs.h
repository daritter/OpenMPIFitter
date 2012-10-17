#ifndef MPIFitter_DspDsmKsPDF_h
#define MPIFitter_DspDsmKsPDF_h

#include <functional>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <boost/foreach.hpp>
//#include "func.h"

#include <Functions.h>
#include "Event.h"
#include <TChain.h>


/** Define all parameters in a namespace to encapsule them with a qualifier.
 *
 * The macro PARAM will define a new parameter with the given name and assign it a unique index as well as
 * fill a map with all known parameter names
 */
namespace PAR {
    PARAM(sig_N);
    PARAM(sig_ratio);
    PARAM(sig_Mbc_mean);
    PARAM(sig_Mbc_sigma);
    PARAM(sig_Mbc_argusC);
    PARAM(sig_dE_ratio);
    PARAM(sig_dE_mean);
    PARAM(sig_dE_meanshift);
    PARAM(sig_dE_sigma);
    PARAM(sig_dE_sigmascale);
    PARAM(sig_dE_cheb1);
    PARAM(bkg_N);
    PARAM(gen_N);
};

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
        CMP_signal      = 1<<0,
        CMP_background  = 1<<1,
        CMP_generic     = 1<<2,
        CMP_continuum   = 1<<3,
        CMP_deltaT      = 1<<4,
        CMP_ALL         = CMP_signal | CMP_background | CMP_generic | CMP_continuum | CMP_deltaT
    };

    DspDsmKsPDF(double lowerMbc, double upperMbc, double lowerdE, double upperdE, const std::vector<std::string> filenames,
            EnabledComponents components=CMP_ALL, int mcInfoRequired=0, int maxPrintOrder = 0):
        components(components), mcInfoRequired(mcInfoRequired),
        maxPrintOrder(maxPrintOrder), nCalls(0), filenames(filenames),
        signal(lowerMbc,upperMbc, lowerdE, upperdE) {}//, background(lowerMbc,upperMbc,lowerdE,upperdE) {}

    /** Return the pdf value for a given paramater set and event */
    double PDF(const DspDsmKsEvent& e, const std::vector<double> &par) const {
        //Set Parameters for signal component
        signal.set(par[PAR::sig_ratio]);
        signal.fcn1.fcnx.set(par[PAR::sig_Mbc_mean], par[PAR::sig_Mbc_sigma]);
        signal.fcn1.fcny.set(par[PAR::sig_dE_ratio],
                par[PAR::sig_dE_mean],  par[PAR::sig_dE_meanshift],
                par[PAR::sig_dE_sigma], par[PAR::sig_dE_sigmascale]);
        signal.fcn2.fcnx.set(e.benergy, par[PAR::sig_Mbc_argusC]);
        signal.fcn2.fcny.set(par[PAR::sig_dE_cheb1]);

        double sig(0), bkg(0), generic(0), continuum(0), deltaT(0);
        if(components & CMP_signal){
            sig = par[PAR::sig_N] * signal(e.Mbc, e.dE);
            //if(components & CMP_signal_Mbc) sig *= signal_Mbc(e.Mbc);
            //if(components & CMP_signal_dE)  sig *= signal_dE(e.dE);
        }
        if(components & CMP_background){
            bkg = par[PAR::bkg_N];
            //FIXME
        }
        if(components & CMP_generic){
            generic = par[PAR::gen_N];
        }
        //FIXME: ...
        return sig + bkg + generic + continuum + deltaT;
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

    /** finalize the event after all processes are collected */
    double finalize(const std::vector<double> &par, double value) const {
        if(components & CMP_signal){
            value -= par[PAR::sig_N];
        }
        if(components & CMP_background){
            value -= par[PAR::bkg_N];
        }
        if(components & CMP_generic){
            value -= par[PAR::gen_N];
        }
        const double logL = -2.0*value;

        //Determine wether to show value
        nCalls++;
        const int order = (nCalls == 0) ? 1 : std::max(std::min((int)std::log10(nCalls), maxPrintOrder), 0);
        const int interval = static_cast<int>(std::pow(10., order));
        if(nCalls % interval == 0){
            std::cout << "call #" << std::setw(5) << std::left << nCalls << ": ";
            std::cout << "-2logL =" << std::setw(18) << std::setprecision(10) << std::scientific << std::right << logL << std::endl;
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
        DspDsmKsEvent event;
        event.setBranches(chain,"bestB");
        for(unsigned int i=start; i<end; ++i){
            if(!chain->GetEntry(i)) {
                std::cerr << "ARRRRRRRRR: There be a problem readin in event " << i << std::endl;
                std::exit(5);
            }
            if(mcInfoRequired>0 && ((event.mcInfo & mcInfoRequired) != mcInfoRequired)) continue;
            data.push_back(event);
        }

        std::cout << "Aye, process " << process << " fully loaded " << data.size()
            << " events and is ready for pillaging" << std::endl;
        delete chain;
    }

    const std::vector<DspDsmKsEvent>& getData() const { return data; }

    protected:

    /** Which components to include into the pdf */
    EnabledComponents components;
    /** Which flag is required from events to be used */
    int mcInfoRequired;
    /** Max order to print out log2L */
    int maxPrintOrder;
    /** Number of calls */
    mutable int nCalls;
    /** vector containing the data values */
    std::vector<DspDsmKsEvent> data;
    /** filename from which the data should be read */
    std::vector<std::string> filenames;

    /** PDF function components */
    mutable Add2DFcn<CompoundFcn2D<Gauss, DoubleGauss>, CompoundFcn2D<Argus, Chebychev1> > signal;
    //mutable Gauss signal_dE;
};

#endif //MPIFitter_DspDsmKsPDF_h
