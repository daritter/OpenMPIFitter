#include <iostream>
#include <fstream>
#include <MPIFitter.h>
#include <Parameters.h>

namespace PAR {
    PARAM(signal_mean);
    PARAM(signal_sigma);
    PARAM(signal_area);
};

struct FitRoutine {
    template<class FCN> int operator()(FCN &fcn){
        params_t params(10);
        Parameters p;
        std::cout << p << std::endl;
        std::ifstream input("params.txt");
        input >> p;
        input.close();
        std::cout << p << std::endl;

        for(int i=0; i<10; ++i){
            params[i] = i;
            fcn(params);
        }
        return 0;
    }
};

int main(int argc, char* argv[]){
    FitRoutine fitter;
    BasePDF pdf;
    MPIFitter core;
    return core.run(fitter, pdf);
}
