import os
import sys

env = Environment(
    ENV=os.environ,
    CXX="mpic++",
    CXXFLAGS=["-O3", "-g", "-Wall"],
    #LINKFLAGS=["-Wl,--as-needed","-Wl,--strip-all"],
    RPATH = "$LIBPATH",
    CPPPATH=["#include"],
    LIBS=[
        "boost_mpi", "boost_serialization", "boost_program_options",
        "Minuit2"
    ],
)
#env.AppendUnique(LINKFLAGS= "-Wl,--strip-all")
try:
    env.ParseConfig("root-config --libs --cflags --ldflags")
except OSError:
    print "Could not find root-config, please make sure ROOT is set up correctly"
    sys.exit(1)

#Path to boost installation if neccessary
for libpath in ("/remote/pcbelle03/ritter/local/",):
    if os.path.exists(libpath):
        env.AppendUnique(LIBPATH=os.path.join(libpath,"lib"))
        env.AppendUnique(CPPPATH=os.path.join(libpath,"include"))

#Disable SCCS and RCS checks
env.SourceCode('.',None)

#Add build dir
VariantDir('build','.',duplicate=0)

#Export environment
Export('env')

#Run SConscript files
SConscript('build/SConscript')
