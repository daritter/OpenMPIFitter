#include <cstdlib>
#include <iostream>
#include <fstream>
#include <MPIFitter.h>
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

void plot_mBCdE(DspDsmKsPDF& pdf, const std::vector<double>& par, TH2D* h_pdf, TH1D* h_bEnergy, int svdVs, const std::string& name= ""){
    Event e;
    e.svdVs = svdVs-1;
    double integral(0);
    //h_bEnergy->GetEffectiveEntries();
    ProgressBar pbar(h_pdf->GetNbinsX()*h_pdf->GetNbinsY()*h_bEnergy->GetNbinsX());
    std::cout << "Plotting mBCdE for '" << name << "': ";
    for(int ix=0; ix<h_pdf->GetNbinsX(); ++ix){
        e.Mbc = h_pdf->GetXaxis()->GetBinCenter(ix+1);
        for(int iy=0; iy<h_pdf->GetNbinsY(); ++iy){
            e.dE  = h_pdf->GetYaxis()->GetBinCenter(iy+1);
            double pdf_value(0);
            double nEvents(0);
            for(int iz=0; iz<h_bEnergy->GetNbinsX(); ++iz){
                std::cout << ++pbar;
                double n = h_bEnergy->GetBinContent(iz+1);
                e.benergy  = h_bEnergy->GetXaxis()->GetBinCenter(iz+1);
                //Beam energy does not exist, no need to evaluate PDF
                if(n<=0) continue; // || e.benergy<e.Mbc) continue;
                nEvents += n;
                pdf_value += n*pdf.PDF(e,par);
            }
            assert(pdf_value==pdf_value);
            //Only Fill the histogram with positive and non-nan values
            if(pdf_value>0) {
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

struct PlotRoutine {
    /** Set some default options */
    PlotRoutine(): parameterIn("params-in.txt"), rootFile("plots"), plotrange_dT(-20,20),
    bins_mBC(60), bins_dE(50), bins_dT(40), sampling_mBC(5), sampling_dE(5), sampling_dT(5), activeComponents(DspDsmKsPDF::CMP_all), no_pdf(false)
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
    bool no_pdf;

    /** Do the plotting */
    template<class FCN> int operator()(FCN &parallel_pdf){
        const Range range_mBC = parallel_pdf.localFCN().getRange_mBC();
        const Range range_dE = parallel_pdf.localFCN().getRange_dE();
        const Range range_dT = parallel_pdf.localFCN().getRange_dT();

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

        std::string names[] = {"signal","misrecon","mixed","charged","continuum"};
        int components[] = {DspDsmKsPDF::CMP_signal, DspDsmKsPDF::CMP_misrecon, DspDsmKsPDF::CMP_mixed, DspDsmKsPDF::CMP_charged, DspDsmKsPDF::CMP_continuum};

        if(activeComponents & DspDsmKsPDF::CMP_deltat){
            std::vector<double> yields[2] = {std::vector<double>(7,0), std::vector<double>(7,0)};
            for(int rbin=0; rbin<7; ++rbin){
                yields[0][rbin] = parallel_pdf.localFCN().get_yield(par, 1, rbin);
                yields[1][rbin] = parallel_pdf.localFCN().get_yield(par, 2, rbin);
            }
            for(int i=0; i<5; ++i){
                int cmp = components[i];
                std::string name = names[i];
                if(!(cmp & activeComponents)) continue;
                parallel_pdf.setOptions(cmp | DspDsmKsPDF::CMP_deltat);

                TH1D *h_dT_fit[56];
                boost::format dt_hist_name("dT_svd%1%_%2%_rbin%3%_fit_%4%_%5%");
                boost::format dt_hist_title("#Deltat, SVD%1%, %5%, rbin %3%, %2%=%4%");
                TH1D* h_dT_fit_all[2][2][2] = {{{0}}};
                int index(0);
                for(int svd=0; svd<2; ++svd){
                    for(int rbin=0; rbin<7; ++rbin){
                        for(int type=0; type<2; ++type){
                            for(int q=0; q<2; ++q){
                                const std::string hname = (dt_hist_name % (svd+1) % (type>0?"qe":"q") % rbin % (q>0?"p":"m") % name).str();
                                const std::string htitle = (dt_hist_title % (svd+1) % (type>0?"qe":"q") % rbin % (q>0?"p":"m") % name).str();
                                h_dT_fit[index++] = new TH1D(hname.c_str(),htitle.c_str(), bins_dT*sampling_dT, plotrange_dT.vmin, plotrange_dT.vmax);
                                if(!h_dT_fit_all[svd][type][q]){
                                    const std::string hname = (dt_hist_name % (svd+1) % (type>0?"qe":"q") % 7 % (q>0?"p":"m") % name).str();
                                    const std::string htitle = (dt_hist_title % (svd+1) % (type>0?"qe":"q") % 7 % (q>0?"p":"m") % name).str();
                                    h_dT_fit_all[svd][type][q] = new TH1D(hname.c_str(),htitle.c_str(), bins_dT*sampling_dT, plotrange_dT.vmin, plotrange_dT.vmax);
                                }
                            }
                        }
                    }
                }

                for(int svd=0; svd<2; ++svd){
                    std::cout << "Plotting dT for " << name << ", SVD" << (svd+1) << ": ";
                    int svdVs = (svd==0?DspDsmKsPDF::PLT_SVD1:DspDsmKsPDF::PLT_SVD2) | DspDsmKsPDF::PLT_DT;
                    std::vector<double> values(3);
                    values[0] = bins_dT*sampling_dT;
                    values[1] = plotrange_dT.vmin;
                    values[2] = plotrange_dT.vmax;
                    std::vector<double> result = parallel_pdf.plotM(svdVs, values, par);
                    index = 0;
                    for(int rbin=0; rbin<7; ++rbin){
                        for(int type=0; type<2; ++type){
                            for(int q=0; q<2; ++q){
                                TH1D* h = h_dT_fit[index+svd*28];
                                for(int i = 0; i<h->GetNbinsX(); ++i){
                                    h->SetBinContent(i+1, result[index*h->GetNbinsX()+i]);
                                }
                                h->Scale(parallel_pdf.localFCN().get_yield(par, svdVs, rbin)/yields[svd][rbin]);
                                h_dT_fit_all[svd][type][q]->Add(h);
                                index++;
                            }
                        }
                    }
                }
            }
        }

        dTCache::print_stats();

        //Parallel is done now, rest is faster in single core. close all other processes and reload all data
        parallel_pdf.close();
        DspDsmKsPDF &local_pdf = parallel_pdf.localFCN();
        if(parallel_pdf.size()>1) local_pdf.load(0,1);
        TH2D *h_MbcdE_data_svd1 = new TH2D("mbcde_svd1_data", "M_{BC}#DeltaE data, SVD1", bins_mBC, range_mBC.vmin, range_mBC.vmax, bins_dE, range_dE.vmin, range_dE.vmax);
        TH2D *h_MbcdE_data_svd2 = new TH2D("mbcde_svd2_data", "M_{BC}#DeltaE data, SVD2", bins_mBC, range_mBC.vmin, range_mBC.vmax, bins_dE, range_dE.vmin, range_dE.vmax);
        TH1D* h_dT_data_q[2][8][2];
        TH1D* h_dT_data_qe[2][8][2];

        for(int svd=0; svd<2; ++svd){
            for(int rbin=0; rbin<8; ++rbin){
                for(int q=0; q<2; ++q){
                    h_dT_data_q[svd][rbin][q] = new TH1D((boost::format("dT_svd%d_q_rbin%d_data_%s") % (svd+1) % rbin % (q>0?"p":"m")).str().c_str(),
                            (boost::format("#Deltat, SVD%d, rbin %d, q=%+d") % (svd+1) % rbin % (2*q-1)).str().c_str(), bins_dT, plotrange_dT.vmin, plotrange_dT.vmax);
                    h_dT_data_qe[svd][rbin][q] = new TH1D((boost::format("dT_svd%d_qe_rbin%d_data_%s") % (svd+1) % rbin % (q>0?"p":"m")).str().c_str(),
                            (boost::format("#Deltat, SVD%d, rbin %d, q#eta=%+d") % (svd+1) % rbin % (2*q-1)).str().c_str(), bins_dT, plotrange_dT.vmin, plotrange_dT.vmax);
                }
            }
        }

        TH1D* h_rbin_data_svd1 = new TH1D("rbin_svd1","rbin fractions, SVD1", 9, -1, 8);
        TH1D* h_rbin_data_svd2 = new TH1D("rbin_svd2","rbin fractions, SVD2", 9, -1, 8);

        TH1D *h_bEnergy_svd1 = new TH1D("svd1_benergy", "Beamenergy, SVD1", 2000, 5.2862, 5.2905);
        TH1D *h_bEnergy_svd2 = new TH1D("svd2_benergy", "Beamenergy, SVD2", 2000, 5.2862, 5.2905);
        TH3D *h_allEvents = new TH3D("mbcdedt_data", "M_{BC}#DeltaE#Deltat data", 200, range_mBC.vmin, range_mBC.vmax, 200, range_dE.vmin, range_dE.vmax, 200, range_dT.vmin, range_dT.vmax);

        for(const Event& e: local_pdf.getData(0)){
            h_bEnergy_svd1->Fill(e.benergy);
            h_MbcdE_data_svd1->Fill(e.Mbc,e.dE);
            h_allEvents->Fill(e.Mbc,e.dE,e.deltaT);
            h_rbin_data_svd1->Fill(e.rbin);
            h_dT_data_q[0][e.rbin][(e.tag_q+1)/2]->Fill(e.deltaT);
            h_dT_data_qe[0][e.rbin][(e.tag_q*e.eta+1)/2]->Fill(e.deltaT);
            h_dT_data_q[0][7][(e.tag_q+1)/2]->Fill(e.deltaT);
            h_dT_data_qe[0][7][(e.tag_q*e.eta+1)/2]->Fill(e.deltaT);
        }
        for(const Event& e: local_pdf.getData(1)){
            h_bEnergy_svd2->Fill(e.benergy);
            h_MbcdE_data_svd2->Fill(e.Mbc,e.dE);
            h_allEvents->Fill(e.Mbc,e.dE,e.deltaT);
            h_rbin_data_svd2->Fill(e.rbin);
            h_dT_data_q[1][e.rbin][(e.tag_q+1)/2]->Fill(e.deltaT);
            h_dT_data_qe[1][e.rbin][(e.tag_q*e.eta+1)/2]->Fill(e.deltaT);
            h_dT_data_q[1][7][(e.tag_q+1)/2]->Fill(e.deltaT);
            h_dT_data_qe[1][7][(e.tag_q*e.eta+1)/2]->Fill(e.deltaT);
        }
        assert(h_bEnergy_svd1->GetBinContent(0)==0 && h_bEnergy_svd1->GetBinContent(h_bEnergy_svd1->GetNbinsX())==0);
        assert(h_bEnergy_svd2->GetBinContent(0)==0 && h_bEnergy_svd2->GetBinContent(h_bEnergy_svd2->GetNbinsX())==0);
        h_rbin_data_svd1->Sumw2();
        h_rbin_data_svd2->Sumw2();
        h_rbin_data_svd1->Scale(1./h_rbin_data_svd1->GetEffectiveEntries());
        h_rbin_data_svd2->Scale(1./h_rbin_data_svd2->GetEffectiveEntries());

        for(int i=0; i<5; ++i){
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

        if(no_pdf) return 0;
        return std::system(("./python/make_plots.py " + rootFile).c_str());
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

    DspDsmKsPDF pdf(1,false);
    PlotRoutine plotter;
    std::string componentList;

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
            ("override", po::value<std::string>(&plotter.overrideParameters)->default_value(plotter.overrideParameters),
             "Aye, give the order to be overrrridn the parrrameters which are given in a comma separated list of name:value pairs")
            ("plot-output,o", po::value<std::string>(&plotter.rootFile)->default_value(plotter.rootFile),
             "Basename to save the plots")
            ("parameter-out,i", po::value<std::string>(&plotter.parameterIn)->default_value(plotter.parameterIn),
             "Thy file to pillage thy parrrametes from")
            ("plot-mindT", po::value<float>(&plotter.plotrange_dT.vmin)->default_value(plotter.plotrange_dT.vmin),
             "The minimal dT value for the plot")
            ("plot-maxdT", po::value<float>(&plotter.plotrange_dT.vmax)->default_value(plotter.plotrange_dT.vmax),
             "The maximal dT value for the plot")
            ("bins-Mbc", po::value<int>(&plotter.bins_mBC)->default_value(plotter.bins_mBC),
             "Number of Bins per axis for the data")
            ("sampling-Mbc", po::value<int>(&plotter.sampling_mBC)->default_value(plotter.sampling_mBC),
             "sampling for the fit")
            ("bins-dE", po::value<int>(&plotter.bins_dE)->default_value(plotter.bins_dE),
             "Number of Bins per axis for the data")
            ("sampling-dE", po::value<int>(&plotter.sampling_dE)->default_value(plotter.sampling_dE),
             "sampling for the fit")
            ("bins-dT", po::value<int>(&plotter.bins_dT)->default_value(plotter.bins_dT),
             "Number of Bins per axis for the data")
            ("sampling-dT", po::value<int>(&plotter.sampling_dT)->default_value(plotter.sampling_dT),
             "sampling for the fit")
            ("no-pdf", po::bool_switch(&plotter.no_pdf),
             "do not create pdfs, just fill the root file")
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
            plotter.activeComponents = DspDsmKsPDF::getComponents(componentList);
        }catch(std::invalid_argument &e){
            std::cout << e.what() << std::endl;
            return 2;
        }
    }
    pdf.setComponents(plotter.activeComponents);

    MPIFitter core;
    return core.run(plotter, pdf);
}
