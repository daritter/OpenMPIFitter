import os

env = Environment(
    CXX="mpic++",
    CXXFLAGS=["-O3", "-g"],
    #LINKFLAGS=["-Wl,--as-needed"],
    RPATH = "$LIBPATH",
    CPPPATH=["#include"],
)
for libpath in ("/remote/pcbelle03/ritter/local/",):
    if os.path.exists(libpath):
        env.AppendUnique(LIBPATH=os.path.join(libpath,"lib"))
        env.AppendUnique(CPPPATH=os.path.join(libpath,"include"))

env.Program("fitter",["src/fitter.cc"],LIBS=["boost_mpi", "boost_serialization"])
