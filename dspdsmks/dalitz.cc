#include <cstdlib>
#include <iostream>
#include <fstream>
#include <Parameters.h>
#include "DspDsmKs.h"
#include "progress.h"

#include <TFile.h>
#include <TH3D.h>
#include <TH2D.h>
#include <TH1D.h>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

namespace po = boost::program_options;

int main(int argc, char* argv[]){

    DspDsmKsPDF pdf(1,false);
    std::string componentList;
    DspDsmKsPDF::EnabledComponents activeComponents = DspDsmKsPDF::CMP_all;

    Range plotrange_mBC("mBC",0,0), plotrange_dE("dE",0,0), plotrange_dT("dT",0,0);
    int bins_splus(80), bins_sminus(80);
    std::string rootFile("dalitz.root"), parameterIn("ddk-out.par");

    try{
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
            ("plot-output,o", po::value<std::string>(&rootFile)->default_value(rootFile),
             "Basename to save the plots")
            ("parameter-out,i", po::value<std::string>(&parameterIn)->default_value(parameterIn),
             "Thy file to pillage thy parrrametes from")
            ("plot-min-Mbc", po::value<float>(&plotrange_mBC.vmin)->default_value(pdf.getRange_mBC().vmin),
             "The minimal Mbc value for the plot")
            ("plot-max-Mbc", po::value<float>(&plotrange_mBC.vmax)->default_value(pdf.getRange_mBC().vmax),
             "The maximal dT value for the plot")
            ("plot-min-dE", po::value<float>(&plotrange_dE.vmin)->default_value(pdf.getRange_dE().vmin),
             "The minimal dE value for the plot")
            ("plot-max-dE", po::value<float>(&plotrange_dE.vmax)->default_value(pdf.getRange_dE().vmax),
             "The maximal dT value for the plot")
            ("plot-min-dT", po::value<float>(&plotrange_dT.vmin)->default_value(pdf.getRange_dT().vmin),
             "The minimal dT value for the plot")
            ("plot-max-dT", po::value<float>(&plotrange_dT.vmax)->default_value(pdf.getRange_dT().vmax),
             "The maximal dT value for the plot")
            ("bins-splus", po::value<int>(&bins_splus)->default_value(bins_splus),
             "Number of Bins per axis for the data")
            ("bins-sminus", po::value<int>(&bins_sminus)->default_value(bins_sminus),
             "Number of Bins per axis for the data")
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
    pdf.load(0,1);

    Parameters params;
    if(!params.load(parameterIn, "")){
        return 2;
    }
    std::vector<double> par = params.getValues();
    std::vector<double> cov{21972.99589, -19397.81411};
    std::vector<double> cov2{-19397.81411, 490464.2685};

    TFile *r_rootFile = new TFile((rootFile+".root").c_str(),"RECREATE");

    //h_dalitz_svd1 = new TH2D( "dalitz_svd1", "sPlot Dalitz, SVD1", bins_splus, 5., 12., bins_sminus, 5., 12.);
    TH2D* h_dalitz_svd2 = new TH2D( "dalitz_svd2", "Normal Dalitz, SVD2", bins_splus, 5., 12., bins_sminus, 5., 12.);
    TH2D* h_pdalitz_svd2 = new TH2D( "pDalitz_svd2", "pPlot Dalitz, SVD2", bins_splus, 5., 12., bins_sminus, 5., 12.);
    TH2D* h_sdalitz_svd2 = new TH2D( "sDalitz_svd2", "sPlot Dalitz, SVD2", bins_splus, 5., 12., bins_sminus, 5., 12.);
    TH2D* h_sdalitz2_svd2 = new TH2D( "sDalitz2_svd2", "sPlot Dalitz s12,s23, SVD2", bins_splus, 0, -1, bins_sminus, 0, -1.);
    TH2D* h_sdalitz3_svd2 = new TH2D( "sDalitz3_svd2", "sPlot Dalitz s23,s23, SVD2", bins_splus, 0, -1, bins_sminus, 0, -1.);

    /*for(const Event& e: pdf.getData(0)){
        if(!(plotrange_mBC(e.Mbc) && plotrange_dE(e.dE) && plotrange_dT(e.dT))) continue;
    }*/
    const double mB0 = 5.27958;
    const double mDs = 2.01029;
    const double mD0 = 1.86486;
    const double mKm = 0.493677;
    const double mKs = 0.497614;
    for(const Event& e: pdf.getData(1)){
        if(!(plotrange_mBC(e.Mbc) && plotrange_dE(e.dE) && plotrange_dT(e.deltaT))) continue;
        const double sWeight = pdf.get_sWeight(e, par, cov[0], cov[1]);
        //const double sW2 = pdf.get_sWeight(e, par, cov2[0], cov2[1]);
        //std::cout << sWeight << " " << sW2 << " " << (sWeight + sW2) << std::endl;
        const double pWeight = pdf.get_cmpFraction(e,par, DspDsmKsPDF::CMP_signal);
        h_dalitz_svd2->Fill(e.m2DspKs, e.m2DsmKs);
        h_pdalitz_svd2->Fill(e.m2DspKs, e.m2DsmKs, pWeight);
        const double s23 = mB0*mB0 + mDs*mDs + mD0*mD0 + mKm*mKm - e.m2DspKs - e.m2DsmKs;
        h_sdalitz_svd2->Fill(e.m2DspKs, e.m2DsmKs, sWeight);
        h_sdalitz2_svd2->Fill(e.m2DspKs, s23, sWeight);
        h_sdalitz3_svd2->Fill(e.m2DsmKs, s23, sWeight);
    }
    r_rootFile->Write();
    r_rootFile->Close();
}
