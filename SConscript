import os
Import('*')

for filename in Glob('src/fitter/*.cc'):
    basename = os.path.basename(filename.abspath)
    executable = "#fitter-" + os.path.splitext(basename)[0]
    env.Program(executable,[filename] + Glob('src/*.cc'))

env.SConscript('veronika/SConscript')
