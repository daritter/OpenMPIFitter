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

#Disable SCCS and RCS checks
env.SourceCode('.',None)

#Add build dir
VariantDir('build','.',duplicate=0)

#Export environment
Export('env')

#Run SConscript files
SConscript('build/SConscript')
