import re
import string
from collections import defaultdict

re_comment = re.compile(r"""
    (?<!\\)              # negative look-behind to not capture \%
    ((?:\\\\)*)          # but allow even number of backslashes in front
    %.*$                 # find comment char
    """,re.X | re.M)
re_macro = re.compile(r"""
    \\                   # start with backslash
    ([^A-Za-z]           # one letter macros
    |[A-Za-z]+)          # normal macros
    """,re.X)

def LineIter(string):
    """Iterate over the lines in a string

    returns lineNr, index of the first char of this line in string and the line
    """

    line = 1
    start = 0
    length = len(string)
    while(start<length):
        try:
            end = string.index("\n",start)+1
        except ValueError:
            end = length+1
        yield line,start,string[start:end]
        start = end
        line +=1

class LatexMacro(object):
    def __init__(self,source,cmd,line,start,pos):
        self.source = source
        self.cmd = cmd
        self.line = line
        self.start = start
        self.pos = pos
        self.reset()

    def reset(self):
        self.charpos = self.start+self.pos[1]

    def __str__(self):
        return "line %4d,%3d: %s" % (self.line,self.pos[0],self.cmd)

    def skip_whitespace(self):
        c = self.charpos
        s = self.source
        while c<len(s) and s[c] in string.whitespace:
            c+=1
        self.charpos = c

    def get_optional(self,delimiter="[]"):
        return self.get_argument(delimiter,True)

    def get_argument(self,delimiter="{}",optional=False):
        try:
            self.skip_whitespace()
            s = self.source
            start = self.charpos
            end = start

            if(s[start] != delimiter[0]):
                if optional: return None
                macro = re_macro.match(s[start:])
                if macro:
                    self.charpos += macro.end()
                    return macro.group()
                else:
                    self.charpos +=1
                    return s[start]

            level=0
            while level>0 or start==end:
                char = s[end]
                if char==delimiter[0]: level+=1
                if char==delimiter[1]: level-=1
                end+=1

            self.charpos = end
            return s[start+1:end-1]
        except IndexError:
            self.charpos = len(self.source)+1
            return None


class LatexEnv(LatexMacro):
    def __init__(self,*args):
        LatexMacro.__init__(self,*args)

    def reset(self):
        LatexMacro.reset(self)
        self.env = self.get_argument()

    def __str__(self):
        return "line %4d,%3d: %s{%s}" % (self.line,self.pos[0],self.cmd,self.env)

class LatexScanner(object):
    def __init__(self,path=['.']):
        self.macro_hooks = defaultdict(list)
        self.begin_hooks = defaultdict(list)
        self.end_hooks = defaultdict(list)
        self.path = path

    def register_macro(self,name,func,**args):
        self.macro_hooks[name].append((func,args))

    def register_begin(self,name,func,**args):
        self.begin_hooks[name].append((func,args))

    def register_end(self,name,func,**args):
        self.end_hooks[name].append((func,args))

    def parse(self,content):
        """
        Parse tex file

        This is a very simple parser which just checks for macros using regular expressions. For each command we check
        if a hook is set to do something with the command. If so, the hook is called with an LatexMacro instance as
        argument
        """

        #remove comments
        data = re_comment.sub(r'\1',content)

        for lineNr,start,line in LineIter(data):
            if line=="": continue

            for cmd in re_macro.finditer(line):
                if cmd.group(1) in ["begin","end"]:
                    cmd = LatexEnv(data, cmd.group(1), lineNr , start, cmd.span(1))
                else:
                    cmd = LatexMacro(data, cmd.group(1), lineNr , start, cmd.span(1))

                def run_hooks(cmd,key,hooks):
                    if key in hooks:
                        for func,args in hooks[key]:
                            cmd.reset()
                            func(self,cmd,**args)

                run_hooks(cmd,cmd.cmd,self.macro_hooks)
                if(cmd.cmd == "begin"): run_hooks(cmd,cmd.env,self.begin_hooks)
                if(cmd.cmd == "end"):   run_hooks(cmd,cmd.env,self.end_hooks)


if __name__ == "__main__":
    ls = LatexScanner()
    def setlength(scanner, c):
        print c, c.get_optional(), c.get_argument()
    ls.register_macro("setlength",setlength)
    ls.register_begin("fmffile",setlength)
    ls.parse(open("test.tex").read())
