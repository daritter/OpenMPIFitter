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
        std::cout << "parameters for part " << part << " out of " << size << ": ";
        for(int i=0; i<params.size(); ++i){
            std::cout << params[i];
            if(i<params.size()-1) std::cout << ", ";
        }
        std::cout << std::endl;
        return 0;
    }

    /** Load the chunk of data to be used by this part of the pdf */
    void load(int argc, char* argv[], int part, int size) {
        this->part = part;
        this->size = size;
    }

    /** Variables used for output, not required */
    int part;
    int size;
};

#endif //DsDsKsFitter_BasePDF_h
