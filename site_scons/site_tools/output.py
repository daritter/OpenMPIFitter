import sys

if sys.stdout.isatty():
    colors = {
        'red'    : '\033[91m',
        'green'  : '\033[92m',
        'yellow' : '\033[93m',
        'blue'   : '\033[94m',
        'purple' : '\033[95m',
        'cyan'   : '\033[96m',
        'end'    : '\033[0m',
    }
else:
    colors = { 'end':'' }

def output(color,cmd):
    color_start = colors.get(color,'end')
    color_end   = colors['end']
    out = '[ %s%s%s ]' % (color_start, cmd.ljust(10), color_end)
    return out

def setOutput(env):
    env.Replace(
      CCCOMSTR        = '${_printCMD("blue","compiling C")} $SOURCES',
      CXXCOMSTR       = '${_printCMD("blue","compiling C++")} $SOURCES',
      SHCCCOMSTR      = '${_printCMD("blue","compiling shared C")} SOURCES',
      SHCXXCOMSTR     = '${_printCMD("blue","compiling shared C++")} $SOURCES',
      F77COMSTR       = '${_prindtCMD("red","compiling F77")} $SOURCES',
      F90COMSTR       = '${_prindtCMD("red","compiling F90")} $SOURCES',
      F95COMSTR       = '${_prindtCMD("red","compiling F95")} $SOURCES',
      FORTRANCOMSTR   = '${_prindtCMD("red","compiling FORTRAN")} $SOURCES',
      SHF77COMSTR     = '${_prindtCMD("red","compiling shared F77")} $SOURCES',
      SHF90COMSTR     = '${_prindtCMD("red","compiling shared F90")} $SOURCES',
      SHF95COMSTR     = '${_prindtCMD("red","compiling shared F95")} $SOURCES',
      SHFORTRANCOMSTR = '${_prindtCMD("red","compiling shared FORTRAN")} $SOURCES',
      ARCOMSTR        = '${_printCMD("yellow","linking static")} $TARGETS',
      SHLINKCOMSTR    = '${_printCMD("yellow","linking shared")} $TARGETS',
      LINKCOMSTR      = '${_printCMD("yellow","linking program")} $TARGETS',
      ROOTCINTCOMSTR  = '${_printCMD("purple","creating root-dict")} $SOURCES',
      INSTALLSTR      = '${_printCMD("green","installing")} $TARGETS',
    )


def generate(env):
    env['_printCMD'] = output
    env.AddMethod(setOutput,'Quiet')

def exists():
    return True
