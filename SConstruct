import os
import sys

kekcc = os.environ.get("HOSTNAME","").find("kek.jp")>=0

tools = ['default', 'output', 'inkscape', 'pdf']
if kekcc:
    tools = ['default','output']

env = Environment(
    tools=tools,
    ENV=os.environ,
    CXX="mpic++",
    CXXFLAGS=["-O3", "-Wall", "-std=c++0x"],
    #LINKFLAGS=["-Wl,--as-needed","-Wl,--strip-all"],
    RPATH = "$LIBPATH",
    CPPPATH=["#include"],
    LIBS=[
        "boost_mpi", "boost_serialization", "boost_program_options", "boost_regex", "boost_random",
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
belle_top_dir = os.environ.get("BELLE_TOP_DIR","/belle/belle/b20090127_0910")
belle_inc_dir = os.path.join(belle_top_dir, "include")
belle_run_dir = os.environ.get("BELLE_RUN_DIR",os.path.join(belle_top_dir,"x86_64-unknown-linux-gnu/opt"))
belle_lib_dir = os.path.join(belle_run_dir, "lib", "so")
belle_cern    = os.path.join(os.environ.get("CERN_ROOT","/belle/cern/2006"), "lib64")

if not os.path.exists(belle_top_dir):
   print "Could not find the belle library"
   sys.exit(1)

env.Append(
    CPPPATH=[belle_inc_dir],
    LIBPATH=[belle_lib_dir, belle_cern],
    LIBS=["tatami", "tables", "belleutil", "packlib", "mathlib"],
)

#Path to boost installation if neccessary
for libpath in ("/remote/pcbelle03/ritter/local/","/home/belle/veronika/local/"):
    if os.path.exists(libpath):
        env.AppendUnique(LIBPATH=os.path.join(libpath,"lib"))
        env.AppendUnique(CPPPATH=os.path.join(libpath,"include"))

#Disable SCCS and RCS checks
env.SourceCode('.',None)

#Add build dir
VariantDir('build','.',duplicate=0)

#Export environment
Export('env','kekcc')

#Run SConscript files
SConscript('build/SConscript')
