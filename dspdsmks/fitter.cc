#include <iostream>
#include <MPIFitter.h>
#include <Parameters.h>
#include <FitRoutine.h>
//#include <PDF/Example.h>
#include "DspDsmKs.h"

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

namespace po = boost::program_options;

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

    /** Read program options using boost::program_options. Could be anything else */
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

    if(!componentList.empty()){
        try{
            activeComponents = DspDsmKsPDF::getComponents(componentList);
        }catch(std::invalid_argument &e){
            std::cout << e.what() << std::endl;
            return 2;
        }
    }
    pdf.setComponents(activeComponents);

    std::string names[] = {"signal","misrecon","mixed","charged"};
    int components[] = {DspDsmKsPDF::CMP_signal,DspDsmKsPDF::CMP_misrecon,DspDsmKsPDF::CMP_mixed,DspDsmKsPDF::CMP_charged};
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

    /** Call the MPI Fitting core and return the result */
    MPIFitter core;
    return core.run(fitter, pdf);
}
