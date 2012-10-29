from scanner import LatexScanner
import sys
import os
import imp
import copy
import subprocess
import re

re_rerun = re.compile(r'''
    (^LaTeX Warning:.*Rerun)|
    (^Package \w+ Warning:.*Rerun)|
    (^LaTeX Warning:.*\n.*Rerun to get citations correct)
''',re.MULTILINE | re.X)

def unique(seq):
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if x not in seen and not seen_add(x)]

def load_module(name,path="packages"):
    try:
        file, path, descr = imp.find_module(name, [os.path.abspath(os.path.join(os.path.dirname(__file__),path))])
        pymodule = imp.load_module(name, file, path, descr)
        return pymodule
    except ImportError:
        return None

class Latex(LatexScanner):
    def __init__(self,target,source,env):
        LatexScanner.__init__(self)

        self.packages = {}
        self.watchfiles = {}
        self.include_only = None
        self.target = target
        self.base = os.path.splitext(target[0].abspath)[0]
        self.source = source[0]
        self.env = env
        self.path = self.add_path(env.get("TEXINPUTS",".:").split(os.path.pathsep))
        self.depends = []

        self.register_macro("documentclass",self.hook_package,type="class")
        self.register_macro("LoadClass",self.hook_package,type="class")
        self.register_macro("usepackage",self.hook_package)
        self.register_macro("RequirePackage",self.hook_package)
        self.register_macro("input",self.hook_input)
        self.register_macro("include",self.hook_input,include=True)

        for suffix in [".aux",".toc",".lof",".lot"]:
            self.register_watchfile(suffix=suffix)

        self.parse(source[0].get_contents())

    def add_path(self,path,prefix=None,suffix=None,keep_old=True):
        new = []
        last = []
        for p in path:
            if p=="":
                last = [""]
                continue

            if prefix is not None and not os.path.isabs(p):
                p = os.path.join(prefix,p)

            if suffix is not None:
                p = os.path.join(p,suffix)

            new.append(os.path.normpath(p))

        if keep_old: new += path
        base = self.env.Dir(".")
        dirs = unique([str(e) for e in base.Rfindalldirs(tuple(new))])
        return tuple(dirs+last)

    def find_file(self,filename,suffixes=[]):
        candidates = [filename] + ["%s.%s" %(filename,suff.lstrip(".")) for suff in suffixes]
        for filename in candidates:
            node = self.env.FindFile(filename,self.path)
            if node: return node

        return None

    def hook_package(self,s,cmd,type="package"):
        options = cmd.get_optional()
        packages = cmd.get_argument()
        for package in packages.split(","):
            self.register_package(package,type=type)

    def hook_input(self,s,cmd,include=False):
        arg = cmd.get_argument()
        #print "input:", arg
        filename = self.find_file(arg,suffixes=[".tex"])
        if filename is None: return
        if include:
            if self.include_only and not arg in self.include_only: return
            self.register_watchfile(os.path.splitex(filename),suffix=".aux")
        self.register_dependency(filename)
        #print "parsing %s" % filename
        self.parse(filename.get_contents())
        #print "finished %s" % filename

    def register_dependency(self,filename):
        self.depends.append(filename)

    def register_watchfile(self,filename=None,suffix=None,action=None,**argk):
        if filename is None:
            filename = self.base

        if suffix is not None:
            filename += '.' + suffix.lstrip('.')

        node = self.env.File(filename)
        self.env.SideEffect(node,self.target[0])
        self.watchfiles[node] = (action,argk)
        #print "Registered watchfile file '%s'" % filename

    def register_package(self,package,type="package"):
        if package in self.packages: return True
        mod = load_module(package)
        self.packages[package] = mod
        if mod is not None:
            mod.register(self)
            return True

        return False



#Building

    def check_watchfiles(self,initial=False):
        if initial:
            self.signatures = {}

        changed = False
        for node,(action,args) in self.watchfiles.items():
            node.clear_memoized_values()
            node.ninfo = node.new_ninfo()
            signature = node.get_csig()
            if signature != self.signatures.get(node,None):
                #if not initial: print "%s changed, running again" % node
                changed = True
                if not initial and action is not None: action(self,node,**args)
            self.signatures[node] = signature

        return changed

    def build(self, texAction, maxretries = 5):
        self.env["ENV"]["TEXINPUTS"] = os.path.pathsep.join(self.path)
        self.check_watchfiles(True)
        log = self.env.File(self.base + ".log")
        for i in range(maxretries):
            ret = texAction(self.target,self.source,self.env)
            if ret!= 0: return ret
            rerun = self.check_watchfiles()
            if re_rerun.search(log.get_contents()):
                rerun = True
                print "Rerunning because log says so"
            if not rerun: break
