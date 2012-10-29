import os
Import('*')

mpiSource = Glob("src/*.cc")
mpiFitter = env.Library("mpifitter", mpiSource)

for filename in Glob('src/fitter/*.cc'):
    basename = os.path.basename(filename.abspath)
    executable = "#fitter-" + os.path.splitext(basename)[0]
    env.Program(executable,[filename] + mpiFitter)

Export(["env", "mpiFitter", "mpiSource"])
env.SConscript('dspdsmks/SConscript')
env.SConscript('note/SConscript')
env.SConscript('python/SConscript')

testEnv = env.Clone()
testEnv.AppendUnique( CPPPATH = ["gtest/include", "gtest"])
gtest = testEnv.Library("gtest", ["gtest/src/gtest-all.cc", "gtest/src/gtest_main.cc"])
testEnv.Program("unittest", Glob('tests/*.cc') + gtest + mpiFitter)
