#include <iostream>
#include <iomanip>
#include <fstream>
#include <MPIFitter.h>
#include <Parameters.h>
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/MnMinos.h"
#include "Minuit2/MnPrint.h"

//#include <FitRoutine.h>
//#include <PDF/Example.h>
#include "DspDsmKs.h"

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

namespace po = boost::program_options;

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
        params.print(false,true);

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

        dTCache::print_stats();

        if(!success){
            std::cout << ANSI_RED;
            std::cout << "ARRRRRRRRR: Minuit is being a harsh mistress and be not converrrrging" << std::endl;
            std::cout << ANSI_END;
        }else{
            std::cout << "Avast, Minuit be converrrrging." << std::endl;
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


/** This is the main function and will be called on all processes with the same
 * arguments.
 *
 * Here we have to setup everything: Parse the program options and set the
 * corresponding variables in the fitting routine and the PDF if any. PDF and
 * Fit function are created here and passed to the MPIFitter for execution.
 * Loading of the acutal data is performed by the MPIFitter and then the
 * FitRoutine is called with a modified PDF which handles the multiprocessing
 */
int main(int argc, char* argv[]){

    DspDsmKsPDF pdf;
    FitRoutine fitter;
    std::string componentList;
    DspDsmKsPDF::EnabledComponents activeComponents = DspDsmKsPDF::CMP_all;
    int svdFlag(0);

    /** Read program options using boost::program_options. Could be anything else */
    try {
        po::options_description desc("Avast, thy options be:");
        desc.add_options()
            ("help,h", "produce this finely crafted help message")
            ("config,c", po::value<std::string>()->default_value("config.ini"),
             "Config file with standard parrrrameters")
            ("input", po::value<std::vector<std::string> >(&pdf.getFiles())->composing(),
             "Root files containing the data")
            ("cmp", po::value<std::string>(&componentList)->default_value(componentList),
             "Comma separated list of components to use for the fit")
            ("min-Mbc", po::value<float>(&pdf.getRange_mBC().vmin)->default_value(pdf.getRange_mBC().vmin),
             "The minimal Mbc value for the fit")
            ("max-Mbc", po::value<float>(&pdf.getRange_mBC().vmax)->default_value(pdf.getRange_mBC().vmax),
             "The maximal Mbc value for the fit")
            ("min-dE", po::value<float>(&pdf.getRange_dE().vmin)->default_value(pdf.getRange_dE().vmin),
             "The minimal dE value for the fit")
            ("max-dE", po::value<float>(&pdf.getRange_dE().vmax)->default_value(pdf.getRange_dE().vmax),
             "The maximal dE value for the fit")
            ("min-dT", po::value<float>(&pdf.getRange_dT().vmin)->default_value(pdf.getRange_dT().vmin),
             "The minimal dT value for the fit")
            ("max-dT", po::value<float>(&pdf.getRange_dT().vmax)->default_value(pdf.getRange_dT().vmax),
             "The maximal dT value for the fit")
            ("bestB", po::value<std::string>(&pdf.getBestB())->default_value(pdf.getBestB()),
             "BestB Selection method to use")
            ("override", po::value<std::string>(&fitter.overrideParameters)->default_value(fitter.overrideParameters),
             "Aye, give the order to be overrrridn the parrrameters which are given in a comma separated list of name:value pairs")
            ("fix", po::value<std::string>(&fitter.fixParameters)->default_value(fitter.fixParameters),
             "Aye, give the order to be fixin the parrrameters which match against this rrrregular expression")
            ("release", po::value<std::string>(&fitter.releaseParameters)->default_value(fitter.releaseParameters),
             "Aye, give the order to be releasin the parrrameters which match against this rrrregular expression")
            ("parameter-in,i", po::value<std::string>(&fitter.parameterIn)->default_value(fitter.parameterIn),
             "Thy file to pillage thy initial parrrameter guesses from")
            ("parameter-out,o", po::value<std::string>(&fitter.parameterOut)->default_value(fitter.parameterOut),
             "The file thy result should be hidden in")
            ("fit-strategy", po::value<int>(&fitter.fitStrategy)->default_value(fitter.fitStrategy),
             "Declare thy preferrrrred strategy to be fittin with")
            ("tries", po::value<int>(&fitter.maxTries)->default_value(fitter.maxTries),
             "Maximum tries to run the fit when Edm>max")
            ("svd", po::value<int>(&svdFlag)->default_value(Component::BOTH),
             "SVD Version to use, 1=svd1, 2=svd2, 3=both")
            ("veto-min-Mbc", po::value<float>(&pdf.getVeto_mBC().vmin)->default_value(pdf.getVeto_mBC().vmin),
             "The minimal Mbc value to veto")
            ("veto-max-Mbc", po::value<float>(&pdf.getVeto_mBC().vmax)->default_value(pdf.getVeto_mBC().vmax),
             "The maximal Mbc value to veto")
            ("veto-min-dE", po::value<float>(&pdf.getVeto_dE().vmin)->default_value(pdf.getVeto_dE().vmin),
             "The minimal dE value to veto")
            ("veto-max-dE", po::value<float>(&pdf.getVeto_dE().vmax)->default_value(pdf.getVeto_dE().vmax),
             "The maximal dE value to veto")
            ("veto", po::bool_switch(&pdf.getVeto()),
             "Apply the veto")
            ;

        po::variables_map vm;
        po::positional_options_description pod;
        pod.add("input", -1);
        po::store(po::command_line_parser(argc, argv).options(desc).positional(pod).run(), vm);
        std::ifstream config(vm["config"].as<std::string>().c_str());
        if(config.is_open()){
            po::store(po::parse_config_file(config, desc, true), vm);
            config.close();
        }
        if (vm.count("help")) {
            std::cout << desc << "\n";
            return 1;
        }
        po::notify(vm);
    }catch(std::exception &e){
        std::cerr<< e.what() << std::endl;
        return 1;
    }

    if(!componentList.empty()){
        try{
            activeComponents = DspDsmKsPDF::getComponents(componentList);
        }catch(std::invalid_argument &e){
            std::cout << e.what() << std::endl;
            return 2;
        }
    }
    pdf.setComponents(activeComponents);
    pdf.setSVD(svdFlag);

    std::string names[] = {"signal", "misrecon", "bbar", "continuum"};
    int components[] = {DspDsmKsPDF::CMP_signal, DspDsmKsPDF::CMP_misrecon, DspDsmKsPDF::CMP_bbar, DspDsmKsPDF::CMP_continuum};
    for(int i=0; i<4; ++i){
        if(!(activeComponents & components[i])){
            if(!fitter.fixParameters.empty()) fitter.fixParameters += "|";
            std::string name = names[i];
            fitter.fixParameters += name+".*|yield_"+name+".*|ratio_"+name+".*";
        }
    }
    if(!(activeComponents & DspDsmKsPDF::CMP_deltat)){
        if(!fitter.fixParameters.empty()) fitter.fixParameters += "|";
        fitter.fixParameters += ".*_dt_.*";
    }
    if(!(activeComponents & ( DspDsmKsPDF::CMP_bbar |  DspDsmKsPDF::CMP_continuum))){
        if(!fitter.fixParameters.empty()) fitter.fixParameters += "|";
        fitter.fixParameters += "bkg.*";
    }
    if(!(svdFlag & Component::SVD1)){
        if(!fitter.fixParameters.empty()) fitter.fixParameters += "|";
        fitter.fixParameters += ".*svd1.*";
    }
    if(!(svdFlag & Component::SVD2)){
        if(!fitter.fixParameters.empty()) fitter.fixParameters += "|";
        fitter.fixParameters += ".*svd2.*";
    }

    /** Call the MPI Fitting core and return the result */
    MPIFitter core;
    return core.run(fitter, pdf);
}
