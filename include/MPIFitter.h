#ifndef DsDsKsFitter_MPIFitter_h
#define DsDsKsFitter_MPIFitter_h

#include <vector>
#include <functional>
#include <boost/mpi.hpp>
#include <boost/serialization/vector.hpp>
#include <BasePDF.h>

typedef std::vector<double> params_t;
BOOST_IS_MPI_DATATYPE(params_t);
BOOST_CLASS_TRACKING(params_t,track_never);
BOOST_CLASS_IMPLEMENTATION(params_t,object_serializable);
BOOST_IS_BITWISE_SERIALIZABLE(params_t);

enum ProcessStatus { PROCESS_CONTINUE, PROCESS_FINISHED };

template<class PDF> class MPIMaster {
    public:
        MPIMaster(boost::mpi::communicator world):world(world), initialized(false) {}
        ~MPIMaster(){
            ProcessStatus status(PROCESS_FINISHED);
            broadcast(world, status, 0);
        }

        void load(int argc, char* argv[]) {
            pdf.load(argc, argv, world.rank(), world.size());
        }

        double operator()(const params_t &params){
            if(!initialized) {
                broadcast(world, boost::mpi::skeleton(const_cast<params_t&>(params)), 0);
                initialized = true;
            }

            ProcessStatus status(PROCESS_CONTINUE);
            broadcast(world, status, 0);
            broadcast(world, boost::mpi::get_content(params), 0);
            double local_result = pdf(params);
            double result(0);
            typename PDF::operator_type op;
            reduce(world, local_result, result, op, 0);
            return result;
        }

    protected:
        boost::mpi::communicator world;
        PDF pdf;
        bool initialized;
};

template<class PDF> class MPIClient {
    public:
        MPIClient(boost::mpi::communicator world):world(world) {}

        void load(int argc, char* argv[]) {
            pdf.load(argc, argv, world.rank(), world.size());
        }

        void run(){
            params_t params;
            ProcessStatus flag;
            broadcast(world, boost::mpi::skeleton(params), 0);
            boost::mpi::content c = boost::mpi::get_content(params);

            typename PDF::operator_type op;
            while(true){
                broadcast(world, flag, 0);
                if(flag == PROCESS_FINISHED) return;
                broadcast(world, c, 0);
                double result = pdf(params);
                reduce(world, result, op, 0);
            }
        }

    protected:
        boost::mpi::communicator world;
        PDF pdf;
};

template<class PDF> class MPIFitter {
    public:
        template<class FitRoutine> int run(int argc, char* argv[], FitRoutine fitter){
            boost::mpi::environment env(argc, argv);
            boost::mpi::communicator world;

            if( world.rank() > 0 ) {
                MPIClient<PDF> client(world);
                client.load(argc, argv);
                client.run();
                return 0;
            }

            std::cout << "Loading Master in a world of size " << world.size() << std::endl;
            MPIMaster<PDF> master(world);
            master.load(argc, argv);

            std::cout << "Calling Fit routine" << std::endl;
            fitter(master);
        }

};

#endif //DsDsKsFitter_MPIFitter_h
