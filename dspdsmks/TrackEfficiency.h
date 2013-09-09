#ifndef MPIFitter_TrackEfficiency_h
#define MPIFitter_TrackEfficiency_h

#include <string>
#include <vector>
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <limits>

/**
 * Class to calculate Tracking efficiency corrections and systematics according
 * to http://belle.kek.jp/secured/wiki/doku.php?id=software:systematics, Belle
 * Note 1176 and information from Jeremy:
 *
 * For fast tracks and pi0 (>200MeV) the efficiency agrees between MC and data
 * but we have 0.35% systematic per track and 4% systematic per pi0. For lower
 * tracks we have additional systematics and efficiency corrections according
 * to belle not 1176 which are uncorrelated, thus can be added in quadrature.
 *
 * This class reads the data files provied by the Belle Note and applies the
 * suggested procedure.
 */
class TrackEfficiency {
    public:
        struct MomentumBin {
            float pmin;
            float pmax;
            float corr;
            float err_uncorr;
            float err_corr;
            float err_syst;
            unsigned int ntracks;

            bool add(double p){
                if(pmin<=p && p<=pmax){
                    ++ntracks;
                    return true;
                }
                return false;
            }
        };

        TrackEfficiency(): ntracks(0), overall(0), eff(1), err(0) {}

        ~TrackEfficiency() {}

        void init(const std::string &filename, bool charged=true){
            overall = charged?trk_syst:pi0_syst;
            err = overall;

            MomentumBin tmp{0};
            std::string line;
            std::ifstream file(filename.c_str());
            if(!file){
                std::cerr << "Could not open file: " << filename << std::endl;
                std::abort();
            }
            while(!file.eof()){
                std::getline(file, line);
                if(line.empty() || line[0] == '%') continue;
                std::stringstream linebuf(line);
                linebuf >> tmp.pmin >> tmp.pmax >> tmp.corr >> tmp.err_uncorr >> tmp.err_corr >> tmp.err_syst;
                //std::cout << tmp.pmin << ", " << tmp.pmax << ", " << tmp.corr << ", " << tmp.err_uncorr << ", " << tmp.err_corr << ", " << tmp.err_syst << tmp.ntracks << std::endl;
                bins.push_back(tmp);
            }
            //Add high momentum bin
            tmp.pmin = 0.0;
            tmp.pmax = std::numeric_limits<double>::infinity();
            tmp.corr = 1.0;
            tmp.err_uncorr = 0;
            tmp.err_corr = 0;
            tmp.err_syst = 0;
            bins.push_back(tmp);
        }

        void add(double plab) {
            for(auto& bin: bins){
                if(bin.add(plab)) break;
            }
            ++ntracks;
        }

        void calculate() {
            if(!ntracks) return;
            eff = 0;
            err = 0;
            for(const auto& bin: bins){
                //std::cout << bin.pmin << ", " << bin.pmax << ", " << bin.corr << ", " << bin.err_uncorr << ", " << bin.err_corr << ", " << bin.err_syst << ", " << bin.ntracks << std::endl;
                //Weighted average of efficiency correction
                eff += bin.corr * bin.ntracks;
                //And weighted average of errors including correlated terms
                double corr = 0;
                for(const auto& bin2: bins){
                    corr += bin.ntracks * bin2.ntracks * bin.err_corr * bin2.err_corr;
                }
                const double n2 = bin.ntracks*bin.ntracks;
                err += n2*bin.err_uncorr*bin.err_uncorr + corr + n2*bin.err_syst*bin.err_syst;
            }
            eff /= ntracks;
            //Add overall tracking systematic and convert from variance to sigma
            err = std::sqrt(err/ntracks/ntracks + eff*eff*overall*overall);
        }

        double get_eff() { return eff; }
        double get_err() { return err; }


    protected:
        static constexpr float trk_syst = 0.0035;
        static constexpr float pi0_syst = 0.04;
        unsigned int ntracks;
        std::vector<MomentumBin> bins;
        float overall;
        double eff;
        double err;
};

#endif
