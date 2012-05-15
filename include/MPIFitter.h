#ifndef MPIFitter_MPIFitter_h
#define MPIFitter_MPIFitter_h

//#define VERBOSE_TIMING

#include <iomanip>
#include <vector>
#include <functional>
#include <boost/mpi.hpp>
#include <boost/foreach.hpp>
#include <Minuit2/FCNBase.h>
#include <unistd.h>

#ifdef VERBOSE_TIMING
#include <sys/time.h>

double getClock(){
    timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return (ts.tv_sec*1e3) + (ts.tv_nsec*1e-6);
}
#endif

typedef std::pair<size_t, double> ParameterChange;

/** Class functioning as master for the fit process.
 * This class will be passed to Minuit as the actual FCN and is responsible to
 * send the parameters to all clients and combine the results
 */
template<class FCN> class MPIMaster: public ROOT::Minuit2::FCNBase {
    public:
        /** Constructor which save a reference to the fcn */
        MPIMaster(boost::mpi::communicator world, FCN &fcn): world(world), fcn(fcn) {}
        /** Destructor */
        virtual ~MPIMaster(){
            //Let all the others now that the are finished by sending a
            //negative number of changed parameters
            int status(-1);
            broadcast(world, status, 0);
        }

        /** FCN evaluation, called by Minuit2 */
        virtual double operator()(const std::vector<double> &params) const {
#ifdef VERBOSE_TIMING
            double start = getClock();
            double split = start;
            double end;
#endif
            //Make a list of changed parameters: On the first call all
            //parameters will be transmitted and only the changed parameters
            //(mostly 1-2) will be transmitted
            changedParameters.clear();
            lastParameters.reserve(params.size());
            for(size_t i=0; i<params.size(); ++i){
                if(i>=lastParameters.size() || params[i]!=lastParameters[i]){
                    changedParameters.push_back(std::make_pair(i,params[i]));
                }
            }
            lastParameters = params;
            //Send the number of changed parameters to all clients
            int changed = changedParameters.size();
            //std::cout << changed << " Parameters have changed" << std::endl;
            broadcast(world, changed, 0);
#ifdef VERBOSE_TIMING
            end = getClock();
            std::cout << "Broadcast of parameters took " << std::fixed
                      << std::setprecision(4) << (end-split) << "ms"
                      << std::endl;
            split = end;
#endif
            //Send the changed parameters to all clients (as binary buffer)
            broadcast(world, (int*) &changedParameters[0],
                    sizeof(ParameterChange)*changed/sizeof(int), 0);

            //Calculate the local result
            double local_result = fcn(params);
#ifdef VERBOSE_TIMING
            end = getClock();
            std::cout << "Local fcn evaluation took " << std::fixed
                << std::setprecision(4) << (end-split) << "ms"
                << std::endl;
            split = end;
#endif
            double result(0);
            //Reduce the result from all clients to one number (usually by adding them)
            typename FCN::operator_type op;
            reduce(world, local_result, result, op, 0);
#ifdef VERBOSE_TIMING
            end = getClock();
            std::cout << "Getting results took " << std::fixed
                << std::setprecision(4) << (end-split) << "ms"
                << std::endl;
            std::cout << "Total time:  " << std::fixed << std::setprecision(4)
                << (end-start) << "ms" << std::endl;
#endif
            //Return the final result to minuit
            return fcn.finalize(params, result);
        }

        /** Return the error_def to Minuit2 */
        virtual double Up() const {
            return FCN::error_def;
        }

    protected:
        /** Communicator to talk to all clients */
        boost::mpi::communicator world;
        /** Array of changed parameters */
        mutable std::vector<ParameterChange> changedParameters;
        /** The parameters used on the last call */
        mutable std::vector<double> lastParameters;
        /** Reference to the FCN */
        FCN &fcn;
};

/** Class representing an MPI client.
 * This class will wait for parameters and evaluate the FCN for those
 * parameters and return the value
 */
template<class FCN> class MPIClient {
    public:
        /** Constructor which save a reference to the fcn */
        MPIClient(boost::mpi::communicator world, FCN &fcn):world(world), fcn(fcn) {}

        /** Run the client.
         * This is basically an endless loop waiting for new parameters and
         * returning the result
         */
        void run(){
            std::vector<double> params;
            std::vector<ParameterChange> changedParameters;
            int changed;

            typename FCN::operator_type op;
            while(true){
                //Recieve the number of changed parameters
                broadcast(world, changed, 0);
                //Negative number of changed parameters means we're finished
                if(changed < 0) return;
                //Otherwise make sure we have enough elements in the list of changed parameters
                changedParameters.resize(changed);
                params.reserve(changed);
                //And then obtain the parameter changes
                broadcast(world, (int*) &changedParameters[0], sizeof(ParameterChange)*changed/sizeof(int), 0);
                //And apply them to the local parameters
                BOOST_FOREACH(ParameterChange &c, changedParameters){
                    if(c.first>=params.size()) params.resize(c.first+1);
                    params[c.first] = c.second;
                }
                //Then call the FCN and return the result to the Master
                double result = fcn(params);
                reduce(world, result, op, 0);
            }
        }

    protected:
        /** Communicator to talk to the master */
        boost::mpi::communicator world;
        /** Reference to the FCN */
        FCN &fcn;
};

/** This class is the actual MPI Fitter.
 * It takes an FCN and a fit routine. The FCN will be wrapped to be usable by
 * minuit2 and then the fit routine will be called with this FCN on the master process and all
 * multi processing will be handled by the wrapped FCN
 */
class MPIFitter {
    public:
        template<class FitRoutine, class FCN> int run(FitRoutine& fitter, FCN& fcn){
            boost::mpi::environment env;
            boost::mpi::communicator world;

            //Load the correct chunk of data
            fcn.load(world.rank(), world.size());
            //If rank>0 we are a client so
            if( world.rank() > 0 ) {
                //Play nice ...
                int nicelevel = nice(20);
                if(nicelevel<19){
                    std::cerr << "Problem nicing process " << world.rank() << std::endl;
                }
                //And start the client
                MPIClient<FCN> client(world, fcn);
                client.run();
                return 0;
            }

            //Otherwise we are the Master
            std::cout << "Loading Master in a world of size " << world.size()
                << " and calling Fit routine" << std::endl;
            //So wrap the FCN and call the fit routine
            MPIMaster<FCN> master(world, fcn);
            return fitter(master);
        }

};

#endif //MPIFitter_MPIFitter_h
