import os
Import('*')

mpiSource = Glob("src/*.cc")
mpiFitter = env.SharedLibrary("mpifitter", mpiSource)

Export(["env", "mpiFitter"])
ddkLib = env.SConscript('dspdsmks/SConscript')
env.SConscript('python/SConscript', exports=['ddkLib'])

if not kekcc:
    env.SConscript('note/SConscript')
    testEnv = env.Clone()
    testEnv.AppendUnique( CPPPATH = ["gtest/include", "gtest"])
    gtest = testEnv.Library("gtest", ["gtest/src/gtest-all.cc", "gtest/src/gtest_main.cc"])
    testEnv.Program("unittest", Glob('tests/*.cc') + gtest + mpiFitter)
