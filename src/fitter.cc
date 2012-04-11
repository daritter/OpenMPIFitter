#include <iostream>
#include <MPIFitter.h>

struct FitRoutine {
    template<class FCN> int operator()(FCN &fcn){
        params_t params(10);

        for(int i=0; i<10; ++i){
            params[i] = i;
            fcn(params);
        }
    }
};


int main(int argc, char* argv[]){
    FitRoutine fitter;
    MPIFitter<BasePDF> core;
    return core.run(argc, argv, fitter);
}
