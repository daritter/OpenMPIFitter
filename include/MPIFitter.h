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

enum MPIFitterCommand {
    EXIT       =  0,
    PARAMETERS =  1,
    OPERATOR   =  2,
    PLOTTER    =  4,
    OPTIONS    =  8
};

struct MPICommand {
    MPICommand(short cmd=0, unsigned short size=0, short flag=0):cmd(cmd), size(size), flag(flag) {}
    short cmd;
    unsigned size;
    unsigned short flag;
};

void broadcast_cmd(const boost::mpi::communicator &world, MPICommand &cmd, int rank){
    broadcast(world, (int*) &cmd, sizeof(cmd)/sizeof(int), rank);
}

template<class T> void broadcast_vec(const boost::mpi::communicator &world, std::vector<T> &parameters, int rank){
    broadcast(world, (int*) &parameters[0], sizeof(T)*parameters.size()/sizeof(int), rank);
}

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
            close();
        }

        void close(){
            //Let all the others now that the are finished by sending a
            //negative number of changed parameters
            MPICommand cmd(EXIT);
            broadcast_cmd(world, cmd, 0);
        }

        void sendParameters(const std::vector<double> &params, int command = PARAMETERS, bool optional=true) const{
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
            if(changedParameters.size()!=0){
                lastParameters = params;
            }else if(optional) return;
            //Send the number of changed parameters to all clients
            MPICommand cmd(command, changedParameters.size());
            //std::cout << changed << " Parameters have changed" << std::endl;
            broadcast_cmd(world, cmd, 0);
            //Send the changed parameters to all clients (as binary buffer)
            broadcast_vec(world, changedParameters, 0);
#ifdef VERBOSE_TIMING
            end = getClock();
            std::cout << "Broadcast of parameters took " << std::fixed
                      << std::setprecision(4) << (end-split) << "ms"
                      << std::endl;
            split = end;
#endif
        }

        /** FCN evaluation, called by Minuit2 */
        virtual double operator()(const std::vector<double> &params) const {
            sendParameters(params, OPERATOR | PARAMETERS, false);
#ifdef VERBOSE_TIMING
            double start = getClock();
            double split = start;
            double end;
#endif
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

        double plot(int flag, std::vector<double> values, const std::vector<double> &params) const {
            sendParameters(params);
            MPICommand cmd(PLOTTER,values.size(),flag);
            broadcast_cmd(world, cmd, 0);
            broadcast_vec(world, values, 0);
            double result = fcn.plot(flag, values, params);
            int nresults(0);
            std::vector<double> results;
            gather(world, result, results, 0);
            result = 0;
            BOOST_FOREACH(double &r, results){
                if(r!=r) continue;
                result += r;
                ++nresults;
            }
            return result / nresults;
        }

        void setOptions(int flag) {
            MPICommand cmd(OPTIONS,0,flag);
            broadcast_cmd(world, cmd, 0);
            fcn.setOptions(flag);
        }


        /** Return the error_def to Minuit2 */
        virtual double Up() const {
            return FCN::error_def;
        }

        FCN& localFCN() { return fcn; }

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
            MPICommand cmd;

            typename FCN::operator_type op;
            world.barrier();
            while(true){
                //Recieve the number of changed parameters
                broadcast_cmd(world, cmd, 0);
                //std::cout << "Recvieved CMD: (" << cmd.cmd << "," << cmd.size << "," << cmd.flag << ")" << std::endl;
                if(cmd.cmd == EXIT) {
                    //std::cout << "Aye, process " << world.rank() << " be exiting" << std::endl;
                    return;
                }
                if(cmd.cmd & PARAMETERS){
                    //Otherwise make sure we have enough elements in the list of changed parameters
                    changedParameters.resize(cmd.size);
                    params.reserve(cmd.size);
                    //And then obtain the parameter changes
                    broadcast_vec(world, changedParameters, 0);
                    //And apply them to the local parameters
                    BOOST_FOREACH(ParameterChange &c, changedParameters){
                        if(c.first>=params.size()) params.resize(c.first+1);
                        params[c.first] = c.second;
                    }
                }
                if(cmd.cmd & OPERATOR){
                    //Then call the FCN and return the result to the Master
                    double result = fcn(params);
                    reduce(world, result, op, 0);
                }
                if(cmd.cmd & PLOTTER){
                    std::vector<double> values(cmd.size,0);
                    broadcast_vec(world, values, 0);
                    double result = fcn.plot(cmd.flag, values, params);
                    //std::cout << "Client Result: " << result << std::endl;
                    gather(world, result, 0);
                }
                if(cmd.cmd & OPTIONS){
                    //std::cout << "Setting Options " << cmd.flag << std::endl;
                    fcn.setOptions(cmd.flag);
                }
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
        template<class FitRoutine, class FCN> int run(FitRoutine& fitter, FCN& fcn, int niceness=20){
            boost::mpi::environment env;
            boost::mpi::communicator world;

            //Load the correct chunk of data
            fcn.load(world.rank(), world.size());
            //If rank>0 we are a client so
            if( world.rank() > 0 ) {
                //Play nice ...
                int nicelevel = nice(niceness);
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
            world.barrier();
            return fitter(master);
        }

};

#endif //MPIFitter_MPIFitter_h
