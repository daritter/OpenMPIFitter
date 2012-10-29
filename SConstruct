import os
import sys

env = Environment(
    tools=['default','output','inkscape','pdf'],
    ENV=os.environ,
    CXX="mpic++",
    CXXFLAGS=["-O3", "-Wall"],
    #LINKFLAGS=["-Wl,--as-needed","-Wl,--strip-all"],
    RPATH = "$LIBPATH",
    CPPPATH=["#include"],
    LIBS=[
        "boost_mpi", "boost_serialization", "boost_program_options", "boost_regex",
        "Minuit2"
    ],
)
env.Quiet()

#env.AppendUnique(LINKFLAGS= "-Wl,--strip-all")
try:
    env.ParseConfig("root-config --libs --cflags --ldflags")
except OSError:
    print "Could not find root-config, please make sure ROOT is set up correctly"
    sys.exit(1)

#Belle library
if not os.path.exists("/belle/belle/b20090127_0910/include/"):
   print "Could not find the belle library"
   sys.exit(1)

env.Append(
    CPPPATH=["/belle/belle/b20090127_0910/include/"],
    LIBPATH=["/belle/belle/b20090127_0910/x86_64-unknown-linux-gnu/opt/lib/so/", "/belle/cern/2006/lib64/"],
    LIBS=["tatami", "tables", "belleutil", "packlib", "mathlib"],
)

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
