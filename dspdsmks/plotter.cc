#include <iostream>
#include <fstream>
#include <Parameters.h>
#include "DspDsmKs.h"

#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>

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

    //FitRoutine fitter;
    std::string parameterIn;
    std::string rootFile;
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
        ("output,o", po::value<std::string>(&rootFile)->default_value(rootFile),
         "Root file to save the plots")
        ("parameter-in,i", po::value<std::string>(&parameterIn)->default_value(parameterIn),
         "Thy file to pillage thy initial parrrameter guesses from")
        ("mcFlag", po::value<int>(&mcInfoRequired)->default_value(mcInfoRequired),
         "Which mc flags to require from the blasted events")
        //("print,p", po::value<int>(&maxPrintOrder)->default_value(maxPrintOrder),
        // "Only print -2logL for each 10^N call")
        //("fix-parameters", po::value<std::string>(&fitter.fixParameters)->default_value(fitter.fixParameters),
        // "Aye, give the order to be fixin the parrrameters which match against this rrrregular expression")
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

    Parameters params;
    std::ifstream input(parameterIn.c_str());
    if(!input){
        std::cerr << "ARRRRRRRRR: Thy parrrrrameter file could not be opened, abandoning ship" << std::endl;
        return 2;
    }
    input >> params;
    input.close();

    std::vector<double> p = params.getValues();

    DspDsmKsPDF pdf(5.20, 5.30, -0.2, 0.2, files,DspDsmKsPDF::CMP_ALL, mcInfoRequired, maxPrintOrder);
    pdf.load(0,1);

    TFile *r_rootFile = new TFile(rootFile.c_str(),"RECREATE");
    TH2D *h_MbcdE_data = new TH2D("mbcde_data","M_{BC}#DeltaE data", 100,5.2,5.3,100,-0.2,0.2);
    TH2D *h_MbcdE_fit = new TH2D("mbcde_fit","M_{BC}#DeltaE fit", 100,5.2,5.3,100,-0.2,0.2);
    TH1D *h_bEnergy = new TH1D("benergy", "Beamenergy", 500, 0,0);
    h_bEnergy->SetBuffer(10000);
    BOOST_FOREACH(const Event& e, pdf.getData()){
        h_bEnergy->Fill(e.benergy);
        h_MbcdE_data->Fill(e.Mbc,e.dE);
    }
    h_bEnergy->BufferEmpty();
    Event e;
    size_t nEvents = pdf.getData().size();
    for(int ix=0; ix<h_MbcdE_fit->GetNbinsX(); ++ix){
        e.Mbc = h_MbcdE_fit->GetXaxis()->GetBinCenter(ix+1);
        for(int iy=0; iy<h_MbcdE_fit->GetNbinsY(); ++iy){
            e.dE  = h_MbcdE_fit->GetYaxis()->GetBinCenter(iy+1);
            double h_pdf(0);
            for(int iz=0; iz<h_bEnergy->GetNbinsX(); ++iz){
                double n = h_bEnergy->GetBinContent(iz+1);
                if(n<=0) continue;
                e.benergy  = h_bEnergy->GetXaxis()->GetBinCenter(iz+1);
                h_pdf += n*pdf.PDF(e,p);
            }
            if(h_pdf>0 && h_pdf==h_pdf) {
                h_MbcdE_fit->Fill(e.Mbc, e.dE, h_pdf/nEvents);
            }
        }
    }
    h_MbcdE_fit->Scale(h_MbcdE_fit->GetXaxis()->GetBinWidth(1)*h_MbcdE_fit->GetYaxis()->GetBinWidth(1));

    h_MbcdE_data->ProjectionX();
    h_MbcdE_data->ProjectionY();
    h_MbcdE_fit->ProjectionX();
    h_MbcdE_fit->ProjectionY();

    r_rootFile->Write();
    r_rootFile->Close();

    return 0;
}
