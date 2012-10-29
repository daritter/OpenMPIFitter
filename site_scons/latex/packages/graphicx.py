def includegraphics(scanner,cmd):
    options = cmd.get_optional()
    graphic = cmd.get_argument()
    filename = scanner.find_file(graphic,["pdf","png","jpg"])
    if filename is None:
        print "ERROR: Could not find graphics file '%s'" % graphic
        filename = graphic
    scanner.register_dependency(filename)

def register(scanner):
    scanner.register_macro("includegraphics",includegraphics)
