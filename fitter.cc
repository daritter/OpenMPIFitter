#define BOOST_MPI_HOMOGENEOUS
#include <boost/mpi.hpp>
#include <iostream>
#include <string>
#include <functional>
#include <cstdlib>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/optional.hpp>
#include <boost/foreach.hpp>
namespace mpi = boost::mpi;

typedef std::vector<double> params_t;
BOOST_CLASS_TRACKING(params_t,track_never);
BOOST_CLASS_IMPLEMENTATION(params_t,object_serializable);
BOOST_IS_BITWISE_SERIALIZABLE(params_t);

/**
 * Master is rank 0.
 * Master broadcasts skeleton of the parameters
 *
 * While not finished {
 * Master broadcasts flag to everybody
 * If flag is task_finished -> clients exit
 * Otherwise:
 * Master broadcasts params
 * Master + Clients reduce PDF
 * reduced to final result on master
 * }
 */

class Client {
    public:
        enum Status { CONTINUE, FINISHED };

        Client(mpi::communicator world):world(world){}
        void load(const std::string &filename);
        void run();
        double evaluate(const params_t &params);
    protected:
        mpi::communicator world;
        std::vector<double> data;
};

class Master {
    public:
        Master(mpi::communicator world):world(world),client(world) {}
        ~Master();
        void load(const std::string &filename) { client.load(filename); }
        void initialize(int nParams);
        double operator()(const params_t &params);
    protected:
        mpi::communicator world;
        Client client;
};

void Client::load(const std::string &filename){
    //FIXME: load relevant rootfile chunk:
    //chunksize =  tree->GetEntriesFast()/world.size();
    //tree->LoadEntry(world.rank()*chunksize);
    //for(int i=0; i<chunksize; ++i) ...
}

void Client::run() {
    params_t params;
    Status flag;
    broadcast(world, mpi::skeleton(params), 0);
    mpi::content c = mpi::get_content(params);

    while(true){
        broadcast(world, flag, 0);
        if(flag == Client::FINISHED) return;
        broadcast(world, c, 0);
        double pdf = evaluate(params);
        reduce(world, pdf, std::multiplies<double>(), 0);
    }
}

double Client::evaluate(const params_t &params){
    //FIXME: need some real pdf here :)
    std::cout << "ev " << world.rank() << ", " << params[0] << std::endl;
    return 0;
}

void Master::initialize(int nParams){
    params_t params(nParams);
    broadcast(world, mpi::skeleton(params), 0);
}

double Master::operator()(const params_t &params){
    Client::Status status(Client::CONTINUE);
    broadcast(world, status, 0);
    broadcast(world, mpi::get_content(params), 0);
    double local_pdf = client.evaluate(params);
    double pdf(0);
    reduce(world, local_pdf, pdf, std::multiplies<double>(), 0);
    return pdf;
}

Master::~Master(){
    Client::Status status(Client::FINISHED);
    broadcast(world, status, 0);
}

int main(int argc, char* argv[]){
    mpi::environment env(argc, argv);
    mpi::communicator world;

    if( world.rank() > 0 ) {
        Client client(world);
        client.run();
        return 0;
    }

    params_t params(1);
    Master master(world);
    master.initialize(1);

    for(int i=0; i<10; ++i){
        params[0] = i;
        master(params);
    }
}


