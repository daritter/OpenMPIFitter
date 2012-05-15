#ifndef MPIFitter_ExamplePDF_h
#define MPIFitter_ExamplePDF_h

#include <functional>
#include <string>
#include <vector>
#include <iostream>
#include "func.h"


/** Define all parameters in a namespace to encapsule them with a qualifier.
 *
 * The macro PARAM will define a new parameter with the given name and assign it a unique index as well as
 * fill a map with all known parameter names
 */
namespace PAR {
    PARAM(sig_N);
    PARAM(sig_mean);
    PARAM(sig_sigma);
    PARAM(bkg_N);
    PARAM(bkg_shape);
};

/** Example PDF function.
 *
 * This is a simple 1D, 2components fit from jeremys tutorial
 */
struct ExamplePDF {

    /** Define how different pieces of the pdf should be reduced to one value
     *
     * useful possibilities include:
     * std::plus<double> -> add the value from all processes
     * std::multiplies<double> -> multiply the value from all processes
     */
    typedef std::plus<double> operator_type;

    /** Define the increase of the objective function value which is eqvivalent to 1 sigma */
    static const double error_def = 1.0;

    /** Return the pdf value for a given paramater set and event */
    double PDF(double e, double x, const std::vector<double> &par) const {
        double sig = par[PAR::sig_N]*Belle::gaussian(x,par[PAR::sig_mean],par[PAR::sig_sigma])/
            Belle::norm_gaussian(x_ll, x_ul, par[PAR::sig_mean],par[PAR::sig_sigma]);
        double bkg = par[PAR::bkg_N]*Belle::argus(x,e,par[PAR::bkg_shape])/Belle::norm_argus(x_ll, x_ul, e,par[PAR::bkg_shape]);
        return sig + bkg;
    }

    /** Return the pdf for all events */
    double operator()(const std::vector<double> &par) const {
        assert( par.size() == 5 );
        double log_pdf(0.0);

        for( unsigned int i=0; i<data_x.size(); i++ ){
            const double pdf = PDF(data_x[i],data_y[i], par);
            if( 0.0 < pdf ) log_pdf += log( pdf ) ;
        }
        return log_pdf;
    }

    /** finalize the event after all processes are collected */
    double finalize(const std::vector<double> &par, const double value) const {
        const double logL = -2.0*(value - par[PAR::sig_N] - par[PAR::bkg_N]);
        std::cout << "-2logL = " << std::setprecision(10) << std::scientific << logL << std::endl;
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
        std::ifstream datafile(filename.c_str());
        double dummy_x(0), dummy_y(0);
        do{
            datafile >> dummy_x >> dummy_y;
            data_x.push_back(dummy_x);
            data_y.push_back(dummy_y);
        }while( !datafile.eof() );

        if(size<=1) return;
        //Now cut the values to be only the segment we really need for this process
        size_t events = (int) std::ceil(1.0 * data_x.size() / size);
        size_t start = process*events;
        size_t end = std::min((process+1)*events,data_x.size());
        std::vector<double> tmp_x;
        tmp_x.reserve(end-start);
        tmp_x.insert(tmp_x.end(), data_x.begin()+start, data_x.begin()+end);
        std::swap(tmp_x,data_x);
        std::vector<double> tmp_y;
        tmp_y.reserve(end-start);
        tmp_y.insert(tmp_y.end(), data_y.begin()+start, data_y.begin()+end);
        std::swap(tmp_y,data_y);
        std::cout << "Process " << process << " of " << size
            << " read in " << data_x.size() << " events"
            << " from " << start
            << " to " << end << std::endl;
    }

    /** Set the parameters given */
    void setParams(const std::string &data, double low, double high){
        filename = data;
        x_ll = low;
        x_ul = high;
    }

    /** vector containing the data x values */
    std::vector<double> data_x;
    /** vector containing the data y values */
    std::vector<double> data_y;
    /** lower limit for x used for normalisation */
    double x_ll;
    /** upper limit for x used for normalisation */
    double x_ul;
    /** filename from which the data should be read */
    std::string filename;
};

#endif //MPIFitter_ExamplePDF_h
