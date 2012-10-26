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

void plot_mBCdE(DspDsmKsPDF& pdf, const std::vector<double>& par, TH2D* h_pdf, TH1D* h_bEnergy, int svdVs, const std::string& name= ""){
    pdf.setSVD((Component::EnabledSVD)svdVs);
    Event e;
    e.svdVs = svdVs-1;
    double integral(0);
    size_t nEvents = h_bEnergy->GetEffectiveEntries();
    ProgressBar pbar(h_pdf->GetNbinsX()*h_pdf->GetNbinsY()*h_bEnergy->GetNbinsX());
    std::cout << "Plotting mBCdE for '" << name << "': ";
    for(int ix=0; ix<h_pdf->GetNbinsX(); ++ix){
        e.Mbc = h_pdf->GetXaxis()->GetBinCenter(ix+1);
        for(int iy=0; iy<h_pdf->GetNbinsY(); ++iy){
            e.dE  = h_pdf->GetYaxis()->GetBinCenter(iy+1);
            double pdf_value(0);
            for(int iz=0; iz<h_bEnergy->GetNbinsX(); ++iz){
                std::cout << ++pbar;
                double n = h_bEnergy->GetBinContent(iz+1);
                if(n<=0) continue;
                e.benergy  = h_bEnergy->GetXaxis()->GetBinCenter(iz+1);
                pdf_value += n*pdf.PDF(e,par);
            }
            if(pdf_value>0 && pdf_value==pdf_value) {
                integral += pdf_value/nEvents
                    * h_pdf->GetXaxis()->GetBinWidth(ix+1)
                    * h_pdf->GetYaxis()->GetBinWidth(iy+1);
                h_pdf->Fill(e.Mbc, e.dE, pdf_value/nEvents);
            }
        }
    }
    double yield = pdf.yield(par);
    h_pdf->Scale(h_pdf->GetXaxis()->GetBinWidth(1)*h_pdf->GetYaxis()->GetBinWidth(1));
    std::cout << "mBCdE Integral for '" << name << "' = " << integral << ", yield = " << yield << ", norm = " << (integral/yield) << std::endl;
}

//template<class FCN> void plotPDF(FCN& pdf, const std::vector<double>& par, TH2D* fit, int svdVs, const std::string& name= ""){
    //Event e;
    //e.svdVs = svdVs;
    //int flag = svdVs | DspDsmKsPDF::PLT_MBCDE;
    //double integral(0);
    //std::vector<double> values(2,0);
    //ProgressBar pbar(fit->GetNbinsX()*fit->GetNbinsY());
    ////size_t nEvents = bEnergy->GetEffectiveEntries();
    //for(int ix=0; ix<fit->GetNbinsX(); ++ix){
        //values[0] = fit->GetXaxis()->GetBinCenter(ix+1);
        //for(int iy=0; iy<fit->GetNbinsY(); ++iy){
            //std::cout << pbar(ix*fit->GetNbinsY()+iy);
            //values[1] = fit->GetYaxis()->GetBinCenter(iy+1);
            //double pdf_value = pdf.plot(flag,values,par);
            //if(pdf_value >0 && pdf_value==pdf_value){
                //integral += pdf_value * fit->GetXaxis()->GetBinWidth(ix+1) *  fit->GetYaxis()->GetBinWidth(ix+1);
                //fit->Fill(values[0], values[1], pdf_value);
            //}
        //}
    //}
    ////double yield = pdf.yield(par);
    //fit->Scale(fit->GetXaxis()->GetBinWidth(1)*fit->GetYaxis()->GetBinWidth(1));
    //std::cout << "PDF Integral for '" << name << "' = " << integral << std::endl; // << ", yield = " << yield << ", norm = " << (integral/yield) << std::endl;
//}

template<class FCN> void plotDT(FCN &pdf, const std::vector<double>& par, TH1D* dtpdf, int flavour, int svdVs, const std::string& name = ""){
    Event e;
    double integral(0);
    int flag = svdVs;
    //if(svdVs == Component::SVD1) flag |= DspDsmKsPDF::PLT_SVD1;
    //if(svdVs == Component::SVD2) flag |= DspDsmKsPDF::PLT_SVD2;
    flag |= flavour>0?DspDsmKsPDF::PLT_DT_P:DspDsmKsPDF::PLT_DT_M;
    std::vector<double> values(1,0);
    ProgressBar pbar(dtpdf->GetNbinsX());
    std::cout << "Plotting dT for '" << name << "': ";
    for(int ix=0; ix<dtpdf->GetNbinsX(); ++ix){
        std::cout << ++pbar;
        values[0] = dtpdf->GetXaxis()->GetBinCenter(ix+1);
        //std::cout << "Plotting dT=" << values[0] << std::endl;
        double pdf_value = pdf.plot(flag,values,par);
        if(pdf_value >0 && pdf_value==pdf_value){
            integral += pdf_value * dtpdf->GetBinWidth(ix+1);
            dtpdf->Fill(values[0], pdf_value);
        }
    }
    //double yield = pdf.localFCN().yield(par);
    dtpdf->Scale(dtpdf->GetBinWidth(1));
    std::cout << "dT Integral for '" << name << "' = " << integral << std::endl;
    //<< ", yield = " << yield << ", norm = " << (integral/yield) << std::endl;
}


struct PlotRoutine {
    /** Set some default options */
    PlotRoutine(): parameterIn("params-in.txt"), rootFile("plots"), plotrange_dT(-70,70),
    nBins(40), sampling(5), nBins_dt(40), sampling_dt(5), activeComponents(DspDsmKsPDF::CMP_all)
    {}

    std::string parameterIn;
    std::string rootFile;
    Range plotrange_dT;
    int nBins;
    int sampling;
    int nBins_dt;
    int sampling_dt;
    DspDsmKsPDF::EnabledComponents activeComponents;

    /** Do the plotting */
    template<class FCN> int operator()(FCN &parallel_pdf){
        const Range range_mBC = parallel_pdf.localFCN().getRange_mBC();
        const Range range_dE = parallel_pdf.localFCN().getRange_dE();
        const Range range_dT = parallel_pdf.localFCN().getRange_dT();

        Parameters params;
        std::ifstream input(parameterIn.c_str());
        if(!input){
            std::cerr << "ARRRRRRRRR: Thy parrrrrameter file could not be opened, abandoning ship" << std::endl;
            return 2;
        }
        input >> params;
        input.close();
        std::vector<double> par = params.getValues();


        TFile *r_rootFile = new TFile((rootFile+".root").c_str(),"RECREATE");


        TH2D *total_MbcdE_fit_svd1 = new TH2D("mbcde_svd1_fit",
                "M_{BC}#DeltaE fit, SVD1", nBins*sampling, range_mBC.vmin, range_mBC.vmax, nBins*sampling, range_dE.vmin, range_dE.vmax);
        TH2D *total_MbcdE_fit_svd2 = new TH2D("mbcde_svd2_fit",
                "M_{BC}#DeltaE fit, SVD2", nBins*sampling, range_mBC.vmin, range_mBC.vmax, nBins*sampling, range_dE.vmin, range_dE.vmax);

        std::string names[] = {"signal","mixed","charged"};
        int components[] = {DspDsmKsPDF::CMP_signal,DspDsmKsPDF::CMP_mixed,DspDsmKsPDF::CMP_charged};
        if(activeComponents & DspDsmKsPDF::CMP_deltat){
        for(int i=0; i<3; ++i){
                int cmp = components[i];
                std::string name = names[i];
                if(!(cmp & activeComponents)) continue;
                if(!name.empty()) name = "_"+name;
                parallel_pdf.setOptions(cmp);

                TH1D *h_dT_fit_svd1_p = new TH1D(("dT_svd1_fit_p" + name).c_str(), "#Deltat fit q=-1, SVD1", nBins_dt*sampling_dt, plotrange_dT.vmin, plotrange_dT.vmax);
                TH1D *h_dT_fit_svd1_m = new TH1D(("dT_svd1_fit_m" + name).c_str(), "#Deltat fit q=+1, SVD1", nBins_dt*sampling_dt, plotrange_dT.vmin, plotrange_dT.vmax);
                TH1D *h_dT_fit_svd2_p = new TH1D(("dT_svd2_fit_p" + name).c_str(), "#Deltat fit q=-1, SVD2", nBins_dt*sampling_dt, plotrange_dT.vmin, plotrange_dT.vmax);
                TH1D *h_dT_fit_svd2_m = new TH1D(("dT_svd2_fit_m" + name).c_str(), "#Deltat fit q=+1, SVD2", nBins_dt*sampling_dt, plotrange_dT.vmin, plotrange_dT.vmax);
                plotDT(parallel_pdf,par,h_dT_fit_svd1_p,+1,DspDsmKsPDF::PLT_SVD1, names[i] + ", SVD1");
                plotDT(parallel_pdf,par,h_dT_fit_svd1_m,-1,DspDsmKsPDF::PLT_SVD1, names[i] + ", SVD1");
                plotDT(parallel_pdf,par,h_dT_fit_svd2_p,+1,DspDsmKsPDF::PLT_SVD2, names[i] + ", SVD2");
                plotDT(parallel_pdf,par,h_dT_fit_svd2_m,-1,DspDsmKsPDF::PLT_SVD2, names[i] + ", SVD2");
            }
        }

        //Parallel is done now, rest is faster in single core. close all other processes and reload all data
        parallel_pdf.close();
        DspDsmKsPDF &local_pdf = parallel_pdf.localFCN();
        local_pdf.load(0,1);
        TH2D *h_MbcdE_data_svd1 = new TH2D("mbcde_svd1_data", "M_{BC}#DeltaE data, SVD1", nBins, range_mBC.vmin, range_mBC.vmax, nBins, range_dE.vmin, range_dE.vmax);
        TH2D *h_MbcdE_data_svd2 = new TH2D("mbcde_svd2_data", "M_{BC}#DeltaE data, SVD2", nBins, range_mBC.vmin, range_mBC.vmax, nBins, range_dE.vmin, range_dE.vmax);
        TH1D *h_dT_data_svd1_p = new TH1D("dT_svd1_data_p", "#Deltat data q=-1, SVD1", nBins_dt, plotrange_dT.vmin, plotrange_dT.vmax);
        TH1D *h_dT_data_svd1_m = new TH1D("dT_svd1_data_m", "#Deltat data q=+1, SVD1", nBins_dt, plotrange_dT.vmin, plotrange_dT.vmax);
        TH1D *h_dT_data_svd2_p = new TH1D("dT_svd2_data_p", "#Deltat data q=-1, SVD2", nBins_dt, plotrange_dT.vmin, plotrange_dT.vmax);
        TH1D *h_dT_data_svd2_m = new TH1D("dT_svd2_data_m", "#Deltat data q=+1, SVD2", nBins_dt, plotrange_dT.vmin, plotrange_dT.vmax);

        TH1D *h_bEnergy_svd1 = new TH1D("svd1_benergy", "Beamenergy, SVD1", 2000, 0,0);
        TH1D *h_bEnergy_svd2 = new TH1D("svd2_benergy", "Beamenergy, SVD2", 2000, 0,0);
        h_bEnergy_svd1->SetBuffer(30000);
        h_bEnergy_svd2->SetBuffer(60000);

        BOOST_FOREACH(const Event& e, local_pdf.getData(0)){
            h_bEnergy_svd1->Fill(e.benergy);
            h_MbcdE_data_svd1->Fill(e.Mbc,e.dE);
            if(e.tag_q>0){
                h_dT_data_svd1_p->Fill(e.deltaT);
            } else {
                h_dT_data_svd1_m->Fill(e.deltaT);
            }
        }
        BOOST_FOREACH(const Event& e, local_pdf.getData(1)){
            h_bEnergy_svd2->Fill(e.benergy);
            h_MbcdE_data_svd2->Fill(e.Mbc,e.dE);
            if(e.tag_q>0){
                h_dT_data_svd2_p->Fill(e.deltaT);
            } else {
                h_dT_data_svd2_m->Fill(e.deltaT);
            }
        }
        h_bEnergy_svd1->BufferEmpty();
        h_bEnergy_svd2->BufferEmpty();

        for(int i=0; i<3; ++i){
            int cmp = components[i];
            std::string name = names[i];
            if(!(cmp & activeComponents)) continue;
            if(!name.empty()) name = "_"+name;

            local_pdf.setOptions(cmp);
            TH2D *h_MbcdE_fit_svd1 = new TH2D(("mbcde_svd1_fit" + name).c_str(),
                    "M_{BC}#DeltaE fit, SVD1", nBins*sampling, range_mBC.vmin, range_mBC.vmax, nBins*sampling, range_dE.vmin, range_dE.vmax);
            TH2D *h_MbcdE_fit_svd2 = new TH2D(("mbcde_svd2_fit" + name).c_str(),
                    "M_{BC}#DeltaE fit, SVD2", nBins*sampling, range_mBC.vmin, range_mBC.vmax, nBins*sampling, range_dE.vmin, range_dE.vmax);

            plot_mBCdE(local_pdf,par,h_MbcdE_fit_svd1, h_bEnergy_svd1, Component::SVD1, names[i] + ", SVD1");
            plot_mBCdE(local_pdf,par,h_MbcdE_fit_svd2, h_bEnergy_svd2, Component::SVD2, names[i] + ", SVD2");

            TH2D* h_MbcdE_fit = (TH2D*) h_MbcdE_fit_svd1->Clone(("mbcde_fit"+name).c_str());
            h_MbcdE_fit->Add(h_MbcdE_fit_svd2);
            total_MbcdE_fit_svd1->Add(h_MbcdE_fit_svd1);
            total_MbcdE_fit_svd2->Add(h_MbcdE_fit_svd2);
        }

        TH2D* total_MbcdE_fit = (TH2D*) total_MbcdE_fit_svd1->Clone("mbcde_fit");
        total_MbcdE_fit->Add(total_MbcdE_fit_svd2);
        TH2D* h_MbcdE_data = (TH2D*) h_MbcdE_data_svd1->Clone("mbcde_data");
        h_MbcdE_data->Add(h_MbcdE_data_svd2);

        r_rootFile->Write();
        r_rootFile->Close();

        return std::system(("./make_plots.py " + rootFile).c_str());
        //return 0;
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

    //FitRoutine fitter;
    PlotRoutine plotter;
    std::vector<std::string> files;
    Range range_mBC(5.24,5.3);
    Range range_dE(-0.1,0.1);
    Range range_dT(-70,70);
    std::string bestB("bestLHsig");
    //DspDsmKsPDF::EnabledComponents activeComponents = DspDsmKsPDF::CMP_all;
    std::vector<std::string> componentList;

    /** Read program options using boost::program_options. Could be anything else */
    po::options_description desc("Avast, thy options be:");
    desc.add_options()
        ("help,h", "produce this finely crafted help message")
        ("config,c", po::value<std::string>()->default_value("config.ini"),
         "Config file with standard parrrrameters")
        ("input", po::value<std::vector<std::string> >(&files)->composing(),
         "Root files containing the data")
        ("plot-output,o", po::value<std::string>(&plotter.rootFile)->default_value(plotter.rootFile),
         "Basename to save the plots")
        ("parameter-out,i", po::value<std::string>(&plotter.parameterIn)->default_value(plotter.parameterIn),
         "Thy file to pillage thy parrrametes from")
        ("minMbc", po::value<float>(&range_mBC.vmin)->default_value(range_mBC.vmin),
         "The minimal Mbc value for the fit")
        ("maxMbc", po::value<float>(&range_mBC.vmax)->default_value(range_mBC.vmax),
         "The maximal Mbc value for the fit")
        ("mindE", po::value<float>(&range_dE.vmin)->default_value(range_dE.vmin),
         "The minimal dE value for the fit")
        ("maxdE", po::value<float>(&range_dE.vmax)->default_value(range_dE.vmax),
         "The maximal dE value for the fit")
        ("mindT", po::value<float>(&range_dT.vmin)->default_value(range_dT.vmin),
         "The minimal dT value for the fit")
        ("maxdT", po::value<float>(&range_dT.vmax)->default_value(range_dT.vmax),
         "The maximal dT value for the fit")
        ("plot-mindT", po::value<float>(&plotter.plotrange_dT.vmin)->default_value(plotter.plotrange_dT.vmin),
         "The minimal dT value for the plot")
        ("plot-maxdT", po::value<float>(&plotter.plotrange_dT.vmax)->default_value(plotter.plotrange_dT.vmax),
         "The maximal dT value for the plot")
        ("bestB", po::value<std::string>(&bestB)->default_value(bestB),
         "BestB Selection method to use")
        ("bins", po::value<int>(&plotter.nBins)->default_value(plotter.nBins),
         "Number of Bins per axis for the data")
        ("sampling", po::value<int>(&plotter.sampling)->default_value(plotter.sampling),
         "sampling for the fit")
        ("bins_dt", po::value<int>(&plotter.nBins_dt)->default_value(plotter.nBins_dt),
         "Number of Bins per axis for the data")
        ("sampling_dt", po::value<int>(&plotter.sampling_dt)->default_value(plotter.sampling_dt),
         "sampling for the fit")
        ("cmp", po::value<std::vector<std::string> >(&componentList)->composing(),
         "Components to use for the fit")
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
        plotter.activeComponents = DspDsmKsPDF::getComponents(componentList);
    }

    DspDsmKsPDF pdf(range_mBC, range_dE, range_dT, files, bestB, (DspDsmKsPDF::EnabledComponents)(plotter.activeComponents ^ DspDsmKsPDF::CMP_deltat), 0);
    MPIFitter core;
    return core.run(plotter, pdf);
}
