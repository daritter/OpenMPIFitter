#ifndef MPIFitter_MPIFitter_h
#define MPIFitter_MPIFitter_h

#include <vector>
#include <functional>
#include <boost/mpi.hpp>
#include <boost/foreach.hpp>
#include <Minuit2/FCNBase.h>
#include <unistd.h>

typedef std::pair<size_t, double> ParameterChange;

template<class PDF> class MPIMaster: public ROOT::Minuit2::FCNBase {
    public:
        MPIMaster(boost::mpi::communicator world, PDF &pdf):world(world), pdf(pdf), initialized(false) {}
        virtual ~MPIMaster(){
            int status(-1);
            broadcast(world, status, 0);
        }

        virtual double operator()(const std::vector<double> &params) const {
            changedParameters.clear();
            lastParameters.reserve(params.size());
            for(size_t i=0; i<params.size(); ++i){
                if(i>=lastParameters.size() || params[i]!=lastParameters[i]){
                    changedParameters.push_back(std::make_pair(i,params[i]));
                }
            }
            lastParameters = params;
            int changed = changedParameters.size();
            broadcast(world, changed, 0);
            broadcast(world, (int*) &changedParameters[0], sizeof(ParameterChange)*changed/sizeof(int), 0);

            double local_result = pdf(params);
            double result(0);
            typename PDF::operator_type op;
            reduce(world, local_result, result, op, 0);
            return pdf.finalize(params, result);
        }

        virtual double Up() const {
            return PDF::error_def;
        }

    protected:
        boost::mpi::communicator world;
        mutable std::vector<ParameterChange> changedParameters;
        mutable std::vector<double> lastParameters;
        PDF &pdf;
        mutable bool initialized;
};

template<class PDF> class MPIClient {
    public:
        MPIClient(boost::mpi::communicator world, PDF &pdf):world(world), pdf(pdf) {}

        void run(){
            std::vector<double> params;
            std::vector<ParameterChange> changedParameters;
            int changed;

            typename PDF::operator_type op;
            while(true){
                broadcast(world, changed, 0);
                if(changed < 0) return;
                changedParameters.resize(changed);
                params.reserve(changed);
                broadcast(world, (int*) &changedParameters[0], sizeof(ParameterChange)*changed/sizeof(int), 0);
                BOOST_FOREACH(ParameterChange &c, changedParameters){
                    if(c.first>=params.size()) params.resize(c.first+1);
                    params[c.first] = c.second;
                }
                double result = pdf(params);
                reduce(world, result, op, 0);
            }
        }

    protected:
        boost::mpi::communicator world;
        PDF &pdf;
};

class MPIFitter {
    public:
        template<class FitRoutine, class PDF> int run(FitRoutine& fitter, PDF& pdf){
            boost::mpi::environment env;
            boost::mpi::communicator world;

            pdf.load(world.rank(), world.size());
            if( world.rank() > 0 ) {
                //Play nice ...
                int nicelevel = nice(20);
                if(nicelevel<19){
                    std::cerr << "Problem nicing process " << world.rank() << std::endl;
                }
                MPIClient<PDF> client(world, pdf);
                client.run();
                return 0;
            }

            std::cout << "Loading Master in a world of size " << world.size() << std::endl;
            MPIMaster<PDF> master(world, pdf);

            std::cout << "Calling Fit routine" << std::endl;
            return fitter(master);
        }

};

#endif //MPIFitter_MPIFitter_h
