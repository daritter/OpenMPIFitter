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
    FitRoutine(): parameterIn("params-in.txt"), parameterOut("params-out.txt"), fixParameters(""), fitStrategy(2), useMinos(false),maxTries(10) {}

    /** Do the fitting */
    template<class FCN> int operator()(FCN &fcn){
        //Loading parameters from file
        Parameters params;
        if(!params.load(parameterIn, overrideParameters, fixParameters, releaseParameters)){
            return 2;
        }
        std::cout << "Aye, this be thy initial parrrrrameters: " << std::endl;
        std::cout << params << std::endl;


        bool success(true);
        //Call Minuit2 to do the actual fit
        if(fitStrategy>0){
            std::ios_base::fmtflags originalFormat = std::cout.flags();
            ROOT::Minuit2::MnUserParameters mnParams = params.getMnParams();
            ROOT::Minuit2::MnMigrad migrad(fcn, mnParams, fitStrategy);
            ROOT::Minuit2::FunctionMinimum min = migrad(50000);
            for(int i=1; i<maxTries; ++i){
                if(!min.IsAboveMaxEdm() || !min.State().IsValid()) break;
                mnParams = min.UserParameters();
                std::cout << "Avast, here goes the fittar for attempt " << (i+1) << std::endl;
                ROOT::Minuit2::MnMigrad migrad(fcn, mnParams, fitStrategy);
                min = migrad(50000);
            }

            std::cout.flags(originalFormat);
            std::cout << min << std::endl;
            std::cout << "Function Minimum: " << std::setprecision(10)
                      << min.Fval() << std::endl;
            if(!min.IsValid()){
                std::cout << "ARRRRRRRRR: Minuit is being a harsh mistress and be not converrrrging" << std::endl;
            }else{
                std::cout << "Avast, Minuit be converrrrging." << std::endl;
            }
            params.update(min.UserParameters());
            success &= min.IsValid();

            /*const ROOT::Minuit2::MnUserCovariance &cov = min.UserCovariance();
            for(unsigned int i=0; i<cov.Nrow(); ++i){
                int external =  min.UserParameters().Trafo().ExtOfInt(i);
                std::cout << "COV row/col " << i << " corresponds to " << external << ": " << min.UserParameters().GetName(external) << std::endl;
            }

            for(unsigned int i=0; i<params.size(); ++i){
                if(params[i].isFixed()) continue;
                std::cout << i << " (" << params[i].name() << ") -> " << min.UserParameters().Trafo().IntOfExt(i) << std::endl;
            }*/
        }

        //Print the non fixed parameters and save all parameters to file
        std::cout << "This be thy final non fixed parrrrameters ye lubber:" << std::endl;
        params.print();
        std::cout << std::endl;

        std::ofstream output(parameterOut.c_str());
        if(!output){
            std::cerr << "ARRRRRRRRR: Thy output parrrrrameter file could not be opened.";
        }else{
            output << params;
            output.close();
        }

        return success?0:1;
    }

    /** Filename to read parameters */
    std::string parameterIn;
    /** Filename to write parameters */
    std::string parameterOut;
    /** Regular expression to determine the parameters we want to fix in
     * addition to the ones in the parameter file */
    std::string fixParameters;
    /** Regular expression to determine the parameters we want to release in
     * addition to the ones in the parameter file */
    std::string releaseParameters;
    /** List of <name>:<value> pairs, possible more than one per string
     * separated by comma for parameter ovverrides to dynamically set a
     * parameter to a different value */
    std::string overrideParameters;
    /** Which strategy to use */
    int fitStrategy;
    /** Wether to call minos after the Fit. Not implemented yet */
    bool useMinos;
    int maxTries;
};

#endif // MPIFitter_FITROUTINE_H
