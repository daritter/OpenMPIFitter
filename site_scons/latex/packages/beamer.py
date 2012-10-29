def register(scanner):
    scanner.register_watchfile(suffix=".head")
    scanner.register_watchfile(suffix=".nav")
    scanner.register_watchfile(suffix=".snm")
    scanner.register_package("hyperref")
