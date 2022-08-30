

def steppable_imports():
    imports = '''from cc3d.cpp.PlayerPython import *\nfrom cc3d import CompuCellSetup
from cc3d.core.PySteppables import *\nimport numpy as np
'''
    return imports


def steppable_declaration(step_name, mitosis=False):
    if mitosis:
        stype = "SteppableBasePy"
    else:
        stype = "MitosisSteppableBase"
    return f"class {step_name}Steppable({stype}):\n"


def mitosis_init(frequency):
    return f'''\n\tdef __init__(self,frequency={frequency}):
\t\tMitosisSteppableBase.__init__(self,frequency)\n'''


def steppable_init(frequency, mitosis=False):
    if not mitosis:
        return mitosis_init(frequency)
    return f'''\n\tdef __init__(self, frequency={frequency}):
\t\tSteppableBasePy.__init__(self,frequency)\n'''


def steppable_start():
    return '''
\tdef start(self):
\t\t"""
\t\tCalled before MCS=0 while building the initial simulation
\t\t"""'''


def generate_steppable(step_name, frequency, mitosis, additional_init=None, additional_start=None):
    imports = steppable_imports()
    declare = steppable_declaration(step_name, mitosis=mitosis)
    init = steppable_init(frequency, mitosis=mitosis)
    if additional_init is not None:
        init = add_to_init(init, additional_init)

    start = steppable_start()
    if additional_start is not None:
        start = add_to_start(start, additional_start)




