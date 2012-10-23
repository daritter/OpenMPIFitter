#include <cstdlib>
#include <iostream>
#include <fstream>
#include <Parameters.h>
#include "DspDsmKs.h"

#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

namespace po = boost::program_options;

void plotPDF(const DspDsmKsPDF& pdf, const std::vector<double>& par, TH2D* fit, TH1D* bEnergy, int svdVs, const std::string& name= ""){
    DspDsmKsEvent e;
    e.svdVs = svdVs;
    double integral(0);
    size_t nEvents = bEnergy->GetEffectiveEntries();
    for(int ix=0; ix<fit->GetNbinsX(); ++ix){
        e.Mbc = fit->GetXaxis()->GetBinCenter(ix+1);
        for(int iy=0; iy<fit->GetNbinsY(); ++iy){
            e.dE  = fit->GetYaxis()->GetBinCenter(iy+1);
            double h_pdf(0);
            for(int iz=0; iz<bEnergy->GetNbinsX(); ++iz){
                double n = bEnergy->GetBinContent(iz+1);
                if(n<=0) continue;
                e.benergy  = bEnergy->GetXaxis()->GetBinCenter(iz+1);
                h_pdf += n*pdf.PDF(e,par);
            }
            if(h_pdf>0 && h_pdf==h_pdf) {
                integral += h_pdf/nEvents
                    * fit->GetXaxis()->GetBinWidth(ix+1)
                    * fit->GetYaxis()->GetBinWidth(iy+1);
                fit->Fill(e.Mbc, e.dE, h_pdf/nEvents);
            }
        }
    }
    double yield = pdf.yield(par);
    fit->Scale(fit->GetXaxis()->GetBinWidth(1)*fit->GetYaxis()->GetBinWidth(1));
    std::cout << "PDF Integral for '" << name << "' = " << integral << ", yield = " << yield << ", norm = " << (integral/yield) << std::endl;
}



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
    std::string rootFile("plots");
    std::vector<std::string> files;
    double lowerMbc(5.2);
    double upperMbc(5.3);
    double lowerdE(-0.2);
    double upperdE( 0.2);
    int nBins(50);
    int oversampling(4);
    DspDsmKsPDF::EnabledComponents activeComponents = DspDsmKsPDF::CMP_ALL;
    std::vector<std::string> componentList;

    //FIXME: components

    /** Read program options using boost::program_options. Could be anything else */
    po::options_description desc("Avast, thy options be:");
    desc.add_options()
        ("help,h", "produce this finely crafted help message")
        ("config,c", po::value<std::string>()->default_value("config.ini"),
         "Config file with standard parrrrameters")
        ("input", po::value<std::vector<std::string> >(&files)->composing(),
         "Root files containing the data")
        ("plot-output,o", po::value<std::string>(&rootFile)->default_value(rootFile),
         "Basename to save the plots")
        ("parameter-out,i", po::value<std::string>(&parameterIn)->default_value(parameterIn),
         "Thy file to pillage thy parrrametes from")
        ("minMbc", po::value<double>(&lowerMbc)->default_value(lowerMbc),
         "The minimal Mbc value for the fit")
        ("maxMbc", po::value<double>(&upperMbc)->default_value(upperMbc),
         "The maximal Mbc value for the fit")
        ("mindE", po::value<double>(&lowerdE)->default_value(lowerdE),
         "The minimal dE value for the fit")
        ("maxdE", po::value<double>(&upperdE)->default_value(upperdE),
         "The maximal dE value for the fit")
        ("bins", po::value<int>(&nBins)->default_value(nBins),
         "Number of Bins per axis for the data")
        ("sampling", po::value<int>(&oversampling)->default_value(oversampling),
         "Oversampling for the fit")
        ("cmp", po::value<std::vector<std::string> >(&componentList)->composing(),
         "Components to use for the fit")
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
        po::store(po::parse_config_file(config, desc, true), vm);
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

    std::vector<double> par = params.getValues();

    if(!componentList.empty()){
        activeComponents = DspDsmKsPDF::CMP_NONE;
        BOOST_FOREACH(std::string &component, componentList){
            boost::to_lower(component);
            boost::trim(component);
            if(component == "signal"){
                activeComponents = (DspDsmKsPDF::EnabledComponents) (activeComponents | DspDsmKsPDF::CMP_signal);
            }
            if(component == "mixed"){
                activeComponents = (DspDsmKsPDF::EnabledComponents) (activeComponents | DspDsmKsPDF::CMP_mixed);
            }
            if(component == "charged"){
                activeComponents = (DspDsmKsPDF::EnabledComponents) (activeComponents | DspDsmKsPDF::CMP_charged);
            }
            if(component == "all"){
                activeComponents = (DspDsmKsPDF::EnabledComponents) (activeComponents | DspDsmKsPDF::CMP_ALL);
            }
        }
    }

    DspDsmKsPDF pdf(lowerMbc, upperMbc, lowerdE, upperdE, files, activeComponents, 0);
    pdf.load(0,1);
    TFile *r_rootFile = new TFile((rootFile+".root").c_str(),"RECREATE");
    TH2D *h_MbcdE_data_svd1 = new TH2D("mbcde_svd1_data", "M_{BC}#DeltaE data, SVD1", nBins, lowerMbc, upperMbc, nBins, lowerdE, upperdE);
    TH2D *h_MbcdE_data_svd2 = new TH2D("mbcde_svd2_data", "M_{BC}#DeltaE data, SVD2", nBins, lowerMbc, upperMbc, nBins, lowerdE, upperdE);
    TH1D *h_bEnergy_svd1 = new TH1D("svd1_benergy", "Beamenergy, SVD1", 500, 0,0);
    TH1D *h_bEnergy_svd2 = new TH1D("svd2_benergy", "Beamenergy, SVD2", 500, 0,0);
    h_bEnergy_svd1->SetBuffer(10000);
    h_bEnergy_svd2->SetBuffer(30000);
    BOOST_FOREACH(const DspDsmKsEvent& e, pdf.getData()){
        if(e.svdVs==0){
            h_bEnergy_svd1->Fill(e.benergy);
            h_MbcdE_data_svd1->Fill(e.Mbc,e.dE);
        }else{
            h_bEnergy_svd2->Fill(e.benergy);
            h_MbcdE_data_svd2->Fill(e.Mbc,e.dE);
        }
    }
    h_bEnergy_svd1->BufferEmpty();
    h_bEnergy_svd2->BufferEmpty();

    TH2D *total_MbcdE_fit_svd1 = new TH2D("mbcde_svd1_fit",
            "M_{BC}#DeltaE fit, SVD1", nBins*oversampling, lowerMbc, upperMbc, nBins*oversampling, lowerdE, upperdE);
    TH2D *total_MbcdE_fit_svd2 = new TH2D("mbcde_svd2_fit",
            "M_{BC}#DeltaE fit, SVD2", nBins*oversampling, lowerMbc, upperMbc, nBins*oversampling, lowerdE, upperdE);

    std::string names[] = {"signal","mixed","charged"};
    int components[] = {DspDsmKsPDF::CMP_signal,DspDsmKsPDF::CMP_mixed,DspDsmKsPDF::CMP_charged};
    for(int i=0; i<3; ++i){
        std::string name = names[i];
        if(!name.empty()) name = "_"+name;
        int cmp = components[i];
        if(!(cmp & activeComponents)) continue;
        pdf.setComponents((DspDsmKsPDF::EnabledComponents) cmp);
        TH2D *h_MbcdE_fit_svd1 = new TH2D(("mbcde_svd1_fit" + name).c_str(),
                "M_{BC}#DeltaE fit, SVD1", nBins*oversampling, lowerMbc, upperMbc, nBins*oversampling, lowerdE, upperdE);
        TH2D *h_MbcdE_fit_svd2 = new TH2D(("mbcde_svd2_fit" + name).c_str(),
                "M_{BC}#DeltaE fit, SVD2", nBins*oversampling, lowerMbc, upperMbc, nBins*oversampling, lowerdE, upperdE);
        pdf.setSVD(Component::SVD1);
        plotPDF(pdf,par,h_MbcdE_fit_svd1,h_bEnergy_svd1,0, names[i] + ", SVD1");
        pdf.setSVD(Component::SVD2);
        plotPDF(pdf,par,h_MbcdE_fit_svd2,h_bEnergy_svd2,1, names[i] + ", SVD2");

        TH2D* h_MbcdE_fit = (TH2D*) h_MbcdE_fit_svd1->Clone(("mbcde_fit"+name).c_str());
        h_MbcdE_fit->Add(h_MbcdE_fit_svd2);
        h_MbcdE_fit->ProjectionX();
        h_MbcdE_fit->ProjectionY();
        h_MbcdE_fit_svd1->ProjectionX();
        h_MbcdE_fit_svd1->ProjectionY();
        h_MbcdE_fit_svd2->ProjectionX();
        h_MbcdE_fit_svd2->ProjectionY();

        total_MbcdE_fit_svd1->Add(h_MbcdE_fit_svd1);
        total_MbcdE_fit_svd2->Add(h_MbcdE_fit_svd2);

    }

    TH2D* h_MbcdE_data = (TH2D*) h_MbcdE_data_svd1->Clone("mbcde_data");
    h_MbcdE_data->Add(h_MbcdE_data_svd2);
    h_MbcdE_data->ProjectionX();
    h_MbcdE_data->ProjectionY();
    h_MbcdE_data_svd1->ProjectionX();
    h_MbcdE_data_svd1->ProjectionY();
    h_MbcdE_data_svd2->ProjectionX();
    h_MbcdE_data_svd2->ProjectionY();

    TH2D* total_MbcdE_fit = (TH2D*) total_MbcdE_fit_svd1->Clone("mbcde_fit");
    total_MbcdE_fit->Add(total_MbcdE_fit_svd2);
    total_MbcdE_fit->ProjectionX();
    total_MbcdE_fit->ProjectionY();
    total_MbcdE_fit_svd1->ProjectionX();
    total_MbcdE_fit_svd1->ProjectionY();
    total_MbcdE_fit_svd2->ProjectionX();
    total_MbcdE_fit_svd2->ProjectionY();

    r_rootFile->Write();
    r_rootFile->Close();

    return std::system(("./make_plots.py " + rootFile).c_str());
}
