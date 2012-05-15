import os
Import('*')

mpiFitter = env.Library("mpifitter",Glob('src/*.cc'))

for filename in Glob('src/fitter/*.cc'):
    basename = os.path.basename(filename.abspath)
    executable = "#fitter-" + os.path.splitext(basename)[0]
    env.Program(executable,[filename] + mpiFitter)

Export(["env","mpiFitter"])
env.SConscript('dspdsmks/SConscript')


testEnv = env.Clone()
testEnv.AppendUnique( CPPPATH = ["gtest/include", "gtest"])
gtest = testEnv.Library("gtest", ["gtest/src/gtest-all.cc", "gtest/src/gtest_main.cc"])
testEnv.Program("unittest", Glob('tests/*.cc') + gtest + mpiFitter)
