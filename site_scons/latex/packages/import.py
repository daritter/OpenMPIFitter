def subimport(scanner,cmd,include=False):
    suffix = cmd.get_argument()
    old_path = scanner.path
    scanner.path = scanner.add_path(old_path,suffix=suffix)
    scanner.hook_input(scanner,cmd,include)
    scanner.path = old_path

def register(scanner):
    scanner.register_macro("subimport",subimport)
    scanner.register_macro("subinputfrom",subimport)
    scanner.register_macro("subincludefrom",subimport,include=True)
