from latex import Latex
import subprocess
import os

from SCons.Action import Action
from SCons.Builder import Builder
import copy

def emitter(target,source,env):
    ls = Latex(target,source,env)
    target[0].latex = ls
    return target,source + ls.depends

def runtex(target,source,env):
    #FIXME: for some reason we miss an environment variable tex needs to find the fonts
    os.environ["TEXINPUTS"] = env["ENV"]["TEXINPUTS"]
    shell_env = copy.copy(env["ENV"])
    for var in "HOME USER USERNAME PWD".split():
        shell_env[var] = os.environ[var]

    tex = subprocess.Popen(["pdflatex","--interaction=nonstopmode","--synctex=1",source[0].abspath],
                           env=shell_env,
                           cwd=os.path.dirname(target[0].abspath),
                           stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
    out = tex.communicate()[0]
    if tex.returncode !=0: print out
    return tex.returncode

runtexAction = Action(runtex,'${_printCMD("red","pdflatex")} $TARGET')
mpostAction =  Action('cd $SOURCE.dir  && mpost $SOURCE.file 2>&1 >/dev/null','${_printCMD("purple","metapost")} $SOURCE')

def build(target,source,env):
    return target[0].latex.build(runtexAction)

buildPdfAction = Action(build,None)#,'${_printCMD("red","pdflatex")} $TARGET')

pdflatex = Builder(
    action=buildPdfAction,
    suffix='.pdf',
    src_suffix='.tex',
    emitter=emitter,
)

def generate(env):
    env.Tool('output')
    env['BUILDERS']['Latex'] = pdflatex
    env['ACTIONS'] = {
        'MetaPost':mpostAction
    }

def exists(env):
    return True
