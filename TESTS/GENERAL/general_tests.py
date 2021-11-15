def install():
    print("In install")


method_name = 'install'  # set by the command line options
possibles = globals().copy()
possibles.update(locals())
method = possibles.get(method_name)
if not method:
    raise NotImplementedError("Method %s not implemented" % method_name)
method()
