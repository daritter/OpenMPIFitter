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
    const double yield = pdf.get_yield(par,svdVs);
    h_pdf->Scale(h_pdf->GetXaxis()->GetBinWidth(1)*h_pdf->GetYaxis()->GetBinWidth(1));
    std::cout << "mBCdE Integral for '" << name << "' = " << integral
        << ", yield = " << yield << ", norm = " << (integral/yield) << std::endl;
}

template<class FCN> void plotDT(FCN &pdf, const std::vector<double>& par, TH1D* dtpdf, int flavour, int flag, const std::string& name = ""){
    Event e;
    double integral(0);
    int svdVs = flag & ( DspDsmKsPDF::PLT_SVD1 | DspDsmKsPDF::PLT_SVD2);
    std::vector<double> values(2,0);
    values[1] = flavour;
    ProgressBar pbar(dtpdf->GetNbinsX());
    std::cout << "Plotting dT for '" << name << "': ";
    for(int ix=0; ix<dtpdf->GetNbinsX(); ++ix){
        std::cout << ++pbar;
        values[0] = dtpdf->GetXaxis()->GetBinCenter(ix+1);
        double pdf_value = pdf.plot(flag,values,par);
        if(pdf_value >0 && pdf_value==pdf_value){
            integral += pdf_value * dtpdf->GetBinWidth(ix+1);
            dtpdf->Fill(values[0], pdf_value);
        }
    }
    const double yield = pdf.localFCN().get_yield(par, svdVs);
    dtpdf->Scale(dtpdf->GetBinWidth(1));
    std::cout << "dT Integral for '" << name << "' = " << integral
        << ", yield = " << yield << ", norm = " << (integral/yield) << std::endl;
}


struct PlotRoutine {
    /** Set some default options */
    PlotRoutine(): parameterIn("params-in.txt"), rootFile("plots"), plotrange_dT(-20,20),
    bins_mBC(60), bins_dE(50), bins_dT(40), sampling_mBC(5), sampling_dE(5), sampling_dT(5), activeComponents(DspDsmKsPDF::CMP_all)
    {}

    std::string parameterIn;
    std::string rootFile;
    std::string overrideParameters;
    Range plotrange_dT;
    int bins_mBC;
    int bins_dE;
    int bins_dT;
    int sampling_mBC;
    int sampling_dE;
    int sampling_dT;
    DspDsmKsPDF::EnabledComponents activeComponents;

    /** Do the plotting */
    template<class FCN> int operator()(FCN &parallel_pdf){
        const Range range_mBC = parallel_pdf.localFCN().getRange_mBC();
        const Range range_dE = parallel_pdf.localFCN().getRange_dE();

        Parameters params;
        if(!params.load(parameterIn, overrideParameters)){
            return 2;
        }
        std::vector<double> par = params.getValues();

        TFile *r_rootFile = new TFile((rootFile+".root").c_str(),"RECREATE");
        TH2D *total_MbcdE_fit_svd1 = new TH2D("mbcde_svd1_fit",
                "M_{BC}#DeltaE fit, SVD1", bins_mBC*sampling_mBC, range_mBC.vmin, range_mBC.vmax, bins_dE*sampling_dE, range_dE.vmin, range_dE.vmax);
        TH2D *total_MbcdE_fit_svd2 = new TH2D("mbcde_svd2_fit",
                "M_{BC}#DeltaE fit, SVD2", bins_mBC*sampling_mBC, range_mBC.vmin, range_mBC.vmax, bins_dE*sampling_dE, range_dE.vmin, range_dE.vmax);

        std::string names[] = {"signal","misrecon","mixed","charged"};
        int components[] = {DspDsmKsPDF::CMP_signal,DspDsmKsPDF::CMP_misrecon,DspDsmKsPDF::CMP_mixed,DspDsmKsPDF::CMP_charged};
        if(activeComponents & DspDsmKsPDF::CMP_deltat){
            for(int i=0; i<4; ++i){
                int cmp = components[i];
                std::string name = names[i];
                if(!(cmp & activeComponents)) continue;
                parallel_pdf.setOptions(cmp);

                for(int i=0;i<2;++i){
                    int flag = i==0?DspDsmKsPDF::PLT_DT_Q:DspDsmKsPDF::PLT_DT_QE;
                    std::string qname = i==0?"q":"qe";
                    for(int svd=0; svd<2; ++svd){
                        int svdVs = svd==0?DspDsmKsPDF::PLT_SVD1:DspDsmKsPDF::PLT_SVD2;
                        for(int flavour=-1; flavour<2; flavour+=2){
                            char fname = flavour<0?'m':'p';
                            std::string histname = (boost::format("dT_svd%d_%s_fit_%s_%s") % (svd+1) % qname % fname % name).str();
                            std::string histtitle = (boost::format("#Deltat fit %s=%+d, SVD%d") % qname % flavour % (svd+1)).str();
                            TH1D* h_dT_fit = new TH1D(histname.c_str(), histtitle.c_str(), bins_dT*sampling_dT, plotrange_dT.vmin, plotrange_dT.vmax);
                            plotDT(parallel_pdf,par,h_dT_fit,flavour,svdVs | flag, (boost::format("%s, %s=%+d SVD%d") % name % qname % flavour % (svd+1)).str());
                        }
                    }
                }
            }
        }

        //Parallel is done now, rest is faster in single core. close all other processes and reload all data
        parallel_pdf.close();
        DspDsmKsPDF &local_pdf = parallel_pdf.localFCN();
        if(parallel_pdf.size()>1) local_pdf.load(0,1);
        TH2D *h_MbcdE_data_svd1 = new TH2D("mbcde_svd1_data", "M_{BC}#DeltaE data, SVD1", bins_mBC, range_mBC.vmin, range_mBC.vmax, bins_dE, range_dE.vmin, range_dE.vmax);
        TH2D *h_MbcdE_data_svd2 = new TH2D("mbcde_svd2_data", "M_{BC}#DeltaE data, SVD2", bins_mBC, range_mBC.vmin, range_mBC.vmax, bins_dE, range_dE.vmin, range_dE.vmax);
        TH1D *h_dT_data_svd1_qep = new TH1D("dT_svd1_qe_data_p", "#Deltat data qe=+1, SVD1", bins_dT, plotrange_dT.vmin, plotrange_dT.vmax);
        TH1D *h_dT_data_svd1_qem = new TH1D("dT_svd1_qe_data_m", "#Deltat data qe=-1, SVD1", bins_dT, plotrange_dT.vmin, plotrange_dT.vmax);
        TH1D *h_dT_data_svd2_qep = new TH1D("dT_svd2_qe_data_p", "#Deltat data qe=+1, SVD2", bins_dT, plotrange_dT.vmin, plotrange_dT.vmax);
        TH1D *h_dT_data_svd2_qem = new TH1D("dT_svd2_qe_data_m", "#Deltat data qe=-1, SVD2", bins_dT, plotrange_dT.vmin, plotrange_dT.vmax);
        TH1D *h_dT_data_svd1_qp = new TH1D("dT_svd1_q_data_p", "#Deltat data q=+1, SVD1", bins_dT, plotrange_dT.vmin, plotrange_dT.vmax);
        TH1D *h_dT_data_svd1_qm = new TH1D("dT_svd1_q_data_m", "#Deltat data q=-1, SVD1", bins_dT, plotrange_dT.vmin, plotrange_dT.vmax);
        TH1D *h_dT_data_svd2_qp = new TH1D("dT_svd2_q_data_p", "#Deltat data q=+1, SVD2", bins_dT, plotrange_dT.vmin, plotrange_dT.vmax);
        TH1D *h_dT_data_svd2_qm = new TH1D("dT_svd2_q_data_m", "#Deltat data q=+1, SVD2", bins_dT, plotrange_dT.vmin, plotrange_dT.vmax);

        TH1D *h_bEnergy_svd1 = new TH1D("svd1_benergy", "Beamenergy, SVD1", 2000, 0,0);
        TH1D *h_bEnergy_svd2 = new TH1D("svd2_benergy", "Beamenergy, SVD2", 2000, 0,0);
        h_bEnergy_svd1->SetBuffer(30000);
        h_bEnergy_svd2->SetBuffer(60000);

        BOOST_FOREACH(const Event& e, local_pdf.getData(0)){
            h_bEnergy_svd1->Fill(e.benergy);
            h_MbcdE_data_svd1->Fill(e.Mbc,e.dE);
            ((e.tag_q*e.eta>0)?h_dT_data_svd1_qep:h_dT_data_svd1_qem)->Fill(e.deltaT);
            ((e.tag_q>0)?h_dT_data_svd1_qp:h_dT_data_svd1_qm)->Fill(e.deltaT);
        }
        BOOST_FOREACH(const Event& e, local_pdf.getData(1)){
            h_bEnergy_svd2->Fill(e.benergy);
            h_MbcdE_data_svd2->Fill(e.Mbc,e.dE);
            ((e.tag_q*e.eta>0)?h_dT_data_svd2_qep:h_dT_data_svd2_qem)->Fill(e.deltaT);
            ((e.tag_q>0)?h_dT_data_svd2_qp:h_dT_data_svd2_qm)->Fill(e.deltaT);
        }
        h_bEnergy_svd1->BufferEmpty();
        h_bEnergy_svd2->BufferEmpty();

        for(int i=0; i<4; ++i){
            int cmp = components[i];
            std::string name = names[i];
            if(!(cmp & activeComponents)) continue;
            if(!name.empty()) name = "_"+name;

            local_pdf.setOptions(cmp);
            TH2D *h_MbcdE_fit_svd1 = new TH2D(("mbcde_svd1_fit" + name).c_str(),
                    "M_{BC}#DeltaE fit, SVD1", bins_mBC*sampling_mBC, range_mBC.vmin, range_mBC.vmax, bins_dE*sampling_dE, range_dE.vmin, range_dE.vmax);
            TH2D *h_MbcdE_fit_svd2 = new TH2D(("mbcde_svd2_fit" + name).c_str(),
                    "M_{BC}#DeltaE fit, SVD2", bins_mBC*sampling_mBC, range_mBC.vmin, range_mBC.vmax, bins_dE*sampling_dE, range_dE.vmin, range_dE.vmax);

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

        return std::system(("./python/make_plots.py " + rootFile).c_str());
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
    Range range_dE(-0.15,0.1);
    Range range_dT(-70,70);
    std::string bestB("bestLHsig");
    //DspDsmKsPDF::EnabledComponents activeComponents = DspDsmKsPDF::CMP_all;
    std::string componentList;

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
        ("bins_Mbc", po::value<int>(&plotter.bins_mBC)->default_value(plotter.bins_mBC),
         "Number of Bins per axis for the data")
        ("sampling_Mbc", po::value<int>(&plotter.sampling_mBC)->default_value(plotter.sampling_mBC),
         "sampling for the fit")
        ("bins_dE", po::value<int>(&plotter.bins_dE)->default_value(plotter.bins_dE),
         "Number of Bins per axis for the data")
        ("sampling_dE", po::value<int>(&plotter.sampling_dE)->default_value(plotter.sampling_dE),
         "sampling for the fit")
        ("bins_dT", po::value<int>(&plotter.bins_dT)->default_value(plotter.bins_dT),
         "Number of Bins per axis for the data")
        ("sampling_dT", po::value<int>(&plotter.sampling_dT)->default_value(plotter.sampling_dT),
         "sampling for the fit")
        ("cmp", po::value<std::string>(&componentList)->default_value(componentList),
         "Components to use for the fit")
        ("override", po::value<std::string>(&plotter.overrideParameters)->default_value(plotter.overrideParameters),
         "Aye, give the order to be releasin the parrrameters which match against this rrrregular expression")
        ;

    po::variables_map vm;
    po::positional_options_description pod;
    pod.add("input", -1);
    po::store(po::command_line_parser(argc, argv).options(desc).positional(pod).allow_unregistered().run(), vm);
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
            plotter.activeComponents = DspDsmKsPDF::getComponents(componentList);
        }catch(std::invalid_argument &e){
            std::cout << e.what() << std::endl;
            return 2;
        }
    }

    DspDsmKsPDF pdf(range_mBC, range_dE, range_dT, files, bestB, (DspDsmKsPDF::EnabledComponents)(plotter.activeComponents ^ DspDsmKsPDF::CMP_deltat), 0);
    MPIFitter core;
    return core.run(plotter, pdf);
}
