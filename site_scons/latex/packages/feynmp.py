def run_mpost(scanner,node):
    scanner.env["ACTIONS"]["MetaPost"](None,node,scanner.env)

def fmffile(scanner,cmd):
    filename = cmd.get_argument()
    scanner.register_watchfile(filename,suffix=".mp",action=run_mpost)

def register(scanner):
    scanner.register_begin("fmffile",fmffile)
