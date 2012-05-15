#include <iostream>
#include <MPIFitter.h>
#include <Parameters.h>
#include <FitRoutine.h>
//#include <PDF/Example.h>
#include "DspDsmKs.h"

#include <boost/program_options.hpp>

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

    FitRoutine fitter;
    std::vector<std::string> files;
    int mcInfoRequired(0);
    int maxPrintOrder(3);

    //FIXME: components

    /** Read program options using boost::program_options. Could be anything else */
    po::options_description desc("Avast, thy options be:");
    desc.add_options()
        ("help,h", "produce this finely crafted help message")
        ("config,c", po::value<std::string>()->default_value("config.ini"),
         "Config file with standard parrrrameters")
        ("input", po::value<std::vector<std::string> >(&files)->composing(),
         "Root files containing the data")
        ("parameter-in,i", po::value<std::string>(&fitter.parameterIn)->default_value(fitter.parameterIn),
         "Thy file to pillage thy initial parrrameter guesses from")
        ("parameter-out,o", po::value<std::string>(&fitter.parameterOut)->default_value(fitter.parameterOut),
         "The file thy result should be hidden in")
        ("fit-strategy", po::value<int>(&fitter.fitStrategy)->default_value(fitter.fitStrategy),
         "Declare thy preferrrrred strategy to be fittin with")
        ("mcFlag", po::value<int>(&mcInfoRequired)->default_value(mcInfoRequired),
         "Which mc flags to require from the blasted events")
        ("print,p", po::value<int>(&maxPrintOrder)->default_value(maxPrintOrder),
         "Only print -2logL for each 10^N call")
        ("fix-parameters", po::value<std::string>(&fitter.fixParameters)->default_value(fitter.fixParameters),
         "Aye, give the order to be fixin the parrrameters which match against this rrrregular expression")
        ;

    po::variables_map vm;
    po::positional_options_description pod;
    pod.add("input", -1);
    po::store(po::command_line_parser(argc, argv).options(desc).positional(pod).run(), vm);
    std::ifstream config(vm["config"].as<std::string>().c_str());
    if(config.is_open()){
        po::store(po::parse_config_file(config, desc), vm);
        config.close();
    }
    if (vm.count("help")) {
        std::cout << desc << "\n";
        return 1;
    }
    po::notify(vm);

    DspDsmKsPDF pdf(5.20, 5.30, -0.2, 0.2, files, DspDsmKsPDF::CMP_ALL, mcInfoRequired, maxPrintOrder);

    /** Call the MPI Fitting core and return the result */
    MPIFitter core;
    return core.run(fitter, pdf);
}
