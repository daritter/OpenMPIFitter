env = Environment(
    CXX="mpic++",
    CXXFLAGS=["-O3"]
)

env.Program("fitter",["fitter.cc"],LIBS=["boost_mpi", "boost_serialization"])
