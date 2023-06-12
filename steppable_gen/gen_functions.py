
def _add_to_function(function, extra):
    return function+'\n'+extra+'\n'


def add_to_init(init, additional_init):
    return _add_to_function(init, additional_init)


def add_to_start(start, additional_start):
    return _add_to_function(start, additional_start)

def add_to_step(step, additional_step):
    return _add_to_function(step, additional_step)


def add_to_finish(finish, additional_finish):
    return _add_to_function(finish, additional_finish)

def add_to_on_stop(on_stop, additional_on_stop):
    return _add_to_function(on_stop, additional_on_stop)

def generate_cell_type_loop(ctype, ntabs):
    tab = ''
    for i in range(ntabs):
        tab += '\t'

    return tab + f"for cell in self.cell_list_by_type(self.{ctype.upper()}):\n"


def steppable_imports(user_data="", phenocell_dir=False):
    if not phenocell_dir:
        phenocell_dir = "C:\\PhenoCellPy"
    imports = '''from cc3d.cpp.PlayerPython import *\nfrom cc3d import CompuCellSetup
from cc3d.core.PySteppables import *\nimport numpy as np\n
'''
    phenocell = f'''import sys\n
# IMPORTANT: PhysiCell has a concept of cell phenotype, PhenoCellPy (https://github.com/JulianoGianlupi/PhenoCellPy) 
# has a similar implementation of phenotypes. You should install PhenoCellPy to translate the Phenotypes from PhysiCell.
# Then change the default path used below with your PhenoCellPy's installation directory
sys.path.extend(['{phenocell_dir}'])
global pcp_imp
pcp_imp = False
try:
\timport PhenoCellPy as pcp
\tpcp_imp = True
except:
\tpass\n\n
user_data={user_data}\n\n
'''
    return imports+phenocell


def steppable_declaration(step_name, mitosis=False):
    if not mitosis:
        stype = "SteppableBasePy"
    else:
        stype = "MitosisSteppableBase"
    return f"class {step_name}Steppable({stype}):\n"


def mitosis_init(frequency):
    return f'''\n\tdef __init__(self,frequency={frequency}):
\t\tMitosisSteppableBase.__init__(self,frequency)\n'''


def steppable_init(frequency, mitosis=False):
    if mitosis:
        return mitosis_init(frequency)
    return f'''\n\tdef __init__(self, frequency={frequency}):
\t\tSteppableBasePy.__init__(self,frequency)\n'''


def steppable_start():
    return '''
\tdef start(self):
\t\t"""
\t\tCalled before MCS=0 while building the initial simulation
\t\t"""
\t\tself.pixel_to_space = float(self.get_xml_element('pixel_to_space').cdata)  # pixel/[unit], see xml for units
\t\tself.mcs_to_time = float(self.get_xml_element('mcs_to_time').cdata)  # MCS/[unit], see xml for units'''


def steppable_step():
    step = '''
\tdef step(self, mcs):
\t\t"""
\t\tCalled every frequency MCS while executing the simulation

\t\t:param mcs: current Monte Carlo step
\t\t"""\n
'''
    return step


def steppable_finish():
    finish = '''
\tdef finish(self):
\t\t"""
\t\tCalled after the last MCS to wrap up the simulation. Good place to close files and do post-processing
\t\t"""\n
'''
    return finish


def steppable_on_stop():
    stop = '''
\tdef on_stop(self):
\t\t"""
\t\tCalled if the simulation is stopped before the last MCS
\t\t"""
\t\tself.finish()\n
'''
    return stop

def mitosis_update_attribute():
    update = '''
\tdef update_attributes(self):
\t\tself.parent_cell.targetVolume /= 2.0
\t\tself.clone_parent_2_child()
'''
    return update

def generate_steppable(step_name, frequency, mitosis, minimal=False, already_imports=False, additional_init=None,
                       additional_start=None, additional_step=None, additional_finish=None, additional_on_stop=None,
                       phenocell_dir=False, user_data=""):
    imports = steppable_imports(user_data=user_data, phenocell_dir=phenocell_dir)
    declare = steppable_declaration(step_name, mitosis=mitosis)
    init = steppable_init(frequency, mitosis=mitosis)
    if additional_init is not None:
        init = add_to_init(init, additional_init)

    start = steppable_start()
    if additional_start is not None:
        start = add_to_start(start, additional_start)
    else:
        start += "\n\t\tpass\n"

    step = steppable_step()
    if additional_step is not None:
        step = add_to_step(step, additional_step)
    else:
        step += "\n\t\tpass\n"

    finish = steppable_finish()

    if additional_finish is not None:
        finish = add_to_finish(finish, additional_finish)

    on_stop = steppable_on_stop()
    
    if additional_on_stop is not None:
        on_stop = add_to_on_stop(on_stop, additional_on_stop)

    mitosis_update = mitosis_update_attribute() if mitosis else ''

    if minimal and already_imports:
        return declare+init+start+"\n"
    elif minimal:
        return imports + declare + init + start + "\n"
    elif not already_imports:
        return imports + declare + init + start + step + mitosis_update + finish + on_stop + "\n"
    return declare+init+start+step+ mitosis_update + finish+on_stop+"\n"


if __name__ == "__main__":
    test_step = generate_steppable("test", 1, False)
    test_mit = generate_steppable("mit", 1, True)
    pass
