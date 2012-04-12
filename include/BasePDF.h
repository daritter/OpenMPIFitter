#ifndef DsDsKsFitter_BasePDF_h
#define DsDsKsFitter_BasePDF_h

#include <functional>
#include <string>
#include <vector>
#include <iostream>

/** Base PDF function which always returns 0 and just prints the parameters.
 *
 * Main use is to show the interface required from real PDFs: the operator()
 * and the load method required by the MPIFitter
 */
struct BasePDF {
    /** Define how different pieces of the pdf should be reduced to one value */
    typedef std::multiplies<double> operator_type;

    /** Return the pdf value for a given paramater set */
    double operator()(const std::vector<double> &params) const {
        std::cout << "parameters for process " << process << " out of " << size << ": ";
        for(int i=0; i<params.size(); ++i){
            std::cout << params[i];
            if(i<params.size()-1) std::cout << ", ";
        }
        std::cout << std::endl;
        return 0;
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
        this->process = process;
        this->size = size;
    }

    /** Variables used for output, not required */
    int process;
    int size;
};

#endif //DsDsKsFitter_BasePDF_h
