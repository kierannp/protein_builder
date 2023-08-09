import sys
from .aminos import *
from .all_names import aminos, just_name

def str_to_class(classname):
    return getattr(sys.modules[__name__], classname)

symbol_class = {}
aas = []
for full, triple, single in aminos.keys():
    if ' ' in full:
        f, s = full.split()
        c_name = f[0].upper() + f[1:] + '_' + s[0].upper()+s[1:]
    else:
        c_name = full[0].upper() + full[1:]

    symbol_class[single] = str_to_class(c_name)