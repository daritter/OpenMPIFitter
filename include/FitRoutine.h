#ifndef MPIFitter_FITROUTINE_H
#define MPIFitter_FITROUTINE_H

#include <iostream>
#include <iomanip>
#include <fstream>
#include <Parameters.h>
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnPrint.h"

/**
 * This struct contains the actual fit routine: The code calling Minuit2 with
 * the given PDF. It does not care wether we run in Multiprocessing mode or
 * not, it just gets the pdf, loads the parameters and calls Minuit2.  It will
 * only be executed in the master process
 */
struct FitRoutine {
    /** Set some default options */
    FitRoutine(): parameterIn("params-in.txt"), parameterOut("params-out.txt"), fixParameters(""), fitStrategy(2), useMinos(false) {}

    /** Do the fitting */
    template<class FCN> int operator()(FCN &fcn){
        Parameters params;
        std::ifstream input(parameterIn.c_str());
        if(!input){
            std::cerr << "ARRRRRRRRR: Thy parrrrrameter file could not be opened, abandoning ship" << std::endl;
            return 2;
        }
        input >> params;
        input.close();

        std::cout << "Aye, this be thy initial parrrrrameters: " << std::endl;
        std::cout << params << std::endl;

        std::ios_base::fmtflags originalFormat = std::cout.flags();
        ROOT::Minuit2::MnUserParameters mnParams = params.getMnParams(fixParameters);
        ROOT::Minuit2::MnMigrad migrad(fcn, mnParams, fitStrategy);
        ROOT::Minuit2::FunctionMinimum min = migrad(10000);
        std::cout.flags(originalFormat);

        std::cout << min << std::endl;
        std::cout << "Function Minimum: " << std::setprecision(10)
		  << min.Fval() << std::endl;
        if(!min.IsValid()){
            std::cerr << "ARRRRRRRRR: Minuit is being a harsh mistress and be not converrrrging" << std::endl;
        }else{
            std::cout << "Avast, Minuit be converrrrging." << std::endl;
        }
        params.update(min.UserParameters());
        std::cout << "This be the final parrrrameters ye lubber:" << std::endl;
        std::cout << params << std::endl;

        std::ofstream output(parameterOut.c_str());
        if(!output){
            std::cerr << "ARRRRRRRRR: Thy output parrrrrameter file could not be opened.";
        }else{
            output << params;
            output.close();
        }

        return 0;
    }

    /** Filename to read parameters */
    std::string parameterIn;
    /** Filename to write parameters */
    std::string parameterOut;
    /** Regular expression to determine the parameters we want to fix in addition to the ones in the parameter file */
    std::string fixParameters;
    /** Which strategy to use */
    int fitStrategy;
    /** Wether to call minos after the Fit. Not implemented yet */
    bool useMinos;
};

#endif // MPIFitter_FITROUTINE_H
