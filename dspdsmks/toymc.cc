#include <cstdlib>
#include <iostream>
#include <fstream>
#include <MPIFitter.h>
#include <Parameters.h>
#include "DspDsmKs.h"
#include "progress.h"

#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

namespace po = boost::program_options;

struct ToyMCRoutine {
    /** Set some default options */
    ToyMCRoutine(): parameterIn("params-in.txt"), output("toymc.root"), scan_dT(-4,4), fudge(1.35),
        scansteps_mBC(100), scansteps_dE(100), scansteps_dT(100), seed(0), gsim(false)
    {}

    std::string parameterIn;
    std::string overrideParameters;
    std::string output;
    Range scan_dT;
    double fudge;
    int scansteps_mBC;
    int scansteps_dE;
    int scansteps_dT;
    unsigned int seed;
    bool gsim;
    bool fullgsim;
    std::vector<std::string> templates;

    /** Do the plotting */
    template<class FCN> int operator()(FCN &parallel_pdf){
        const Range range_mBC = parallel_pdf.localFCN().getRange_mBC();
        const Range range_dE = parallel_pdf.localFCN().getRange_dE();

        //If we have templates but no gsim parameter we generate from pdf, thus removing the templates
        if(!gsim) templates.clear();

        Parameters params;
        if(!params.load(parameterIn, overrideParameters)){
            return 2;
        }
        std::vector<double> par = params.getValues();

        //Now we need to get the maximal value of the pdf
        double step_mBC = (range_mBC.vmax-range_mBC.vmin)/scansteps_mBC;
        double step_dE = (range_dE.vmax-range_dE.vmin)/scansteps_dE;
        double step_dT = (scan_dT.vmax-scan_dT.vmin)/scansteps_dT;

        std::vector<double> values(10,0);
        values[0] = 5.289;
        values[1] = range_mBC.vmin;
        values[2] = range_mBC.vmax;
        values[3] = step_mBC;
        values[4] = range_dE.vmin;
        values[5] = range_dE.vmax;
        values[6] = step_dE;
        values[7] = scan_dT.vmin;
        values[8] = scan_dT.vmax;
        values[9] = step_dT;
        double maxVal[2];
        int flag = templates.empty()?DspDsmKsPDF::PLT_MAX:DspDsmKsPDF::PLT_MAXDT;
        maxVal[0] = fudge*parallel_pdf.plot(DspDsmKsPDF::PLT_SVD1 | flag, values, par, OP_MAX);
        maxVal[1] = fudge*parallel_pdf.plot(DspDsmKsPDF::PLT_SVD2 | flag, values, par, OP_MAX);

        //Parallel is done now, rest is only possible in single core. Close all other processes and reload all data
        parallel_pdf.close();
        DspDsmKsPDF &local_pdf = parallel_pdf.localFCN();
        if(parallel_pdf.size()>1) local_pdf.load(0,1);

        TFile* f = new TFile(output.c_str(),"RECREATE");
        TTree* tree = new TTree("B0","B0 Toy MC");
        local_pdf.generateToyMC(tree,par,maxVal,templates,seed,fullgsim);
        tree->Write();
        f->Write();
        f->Close();
        return 0;
    }
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
    ToyMCRoutine toymc;
    DspDsmKsPDF::EnabledComponents activeComponents = DspDsmKsPDF::CMP_all;
    std::string componentList;

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
        ("override", po::value<std::string>(&toymc.overrideParameters)->default_value(toymc.overrideParameters),
         "Aye, give the order to be modyfying parrrameters given with this list")
        ("toymc-output,o", po::value<std::string>(&toymc.output)->default_value(toymc.output),
         "The name of thy rrroot file to store the events generated by this fine programme")
        ("parameter-out,i", po::value<std::string>(&toymc.parameterIn)->default_value(toymc.parameterIn),
         "Thy file to pillage thy parrrametes from")
        ("scan-min-dT", po::value<float>(&toymc.scan_dT.vmin)->default_value(toymc.scan_dT.vmin),
         "Specify thy minimal dT value to scan for thy maximal pdf value")
        ("scan-max-dT", po::value<float>(&toymc.scan_dT.vmax)->default_value(toymc.scan_dT.vmax),
         "Specify thy maximal dT value to scan for thy maximal pdf value")
        ("steps-Mbc", po::value<int>(&toymc.scansteps_mBC)->default_value(toymc.scansteps_mBC),
         "Name the number of steps to be used in scanning thy Mbc range")
        ("steps-dE", po::value<int>(&toymc.scansteps_dE)->default_value(toymc.scansteps_dE),
         "Name the number of steps to be used in scanning thy dE range")
        ("steps-dT", po::value<int>(&toymc.scansteps_dT)->default_value(toymc.scansteps_dT),
         "Name the number of steps to be used in scanning thy dT range")
        ("fudge", po::value<double>(&toymc.fudge)->default_value(toymc.fudge),
         "Fudgefactor we'd be applying to the maximal pdf value to be on the safe side")
        ("seed", po::value<unsigned int>(&toymc.seed)->default_value(toymc.seed),
         "Seed to use for the number generator, 0=initialize from /dev/urandom")
        ("template", po::value<std::vector<std::string> >(&toymc.templates)->composing(),
         "Data template to draw dE/Mbc from")
        ("gsim", po::bool_switch(&toymc.gsim),
         "Wether to generate from PDF or from gsim. If no templates are given we always generate from pdf")
        ("fullgsim", po::bool_switch(&toymc.fullgsim),
         "Wether to also take dT values from gsim (not applicable for signal)")
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

    MPIFitter core;
    return core.run(toymc, pdf);
}
