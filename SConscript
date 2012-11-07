import os
Import('*')

mpiSource = Glob("src/*.cc")
mpiFitter = env.Library("mpifitter", mpiSource)

Export(["env", "mpiFitter", "mpiSource"])
env.SConscript('dspdsmks/SConscript')
env.SConscript('python/SConscript')

if not kekcc:
    env.SConscript('note/SConscript')
    testEnv = env.Clone()
    testEnv.AppendUnique( CPPPATH = ["gtest/include", "gtest"])
    gtest = testEnv.Library("gtest", ["gtest/src/gtest-all.cc", "gtest/src/gtest_main.cc"])
    testEnv.Program("unittest", Glob('tests/*.cc') + gtest + mpiFitter)
