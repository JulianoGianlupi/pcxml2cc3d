from .generate_constraint_step import get_dicts_for_type

try:
    from .gen_functions import generate_steppable, steppable_imports
except:
    from gen_functions import generate_steppable, steppable_imports  # why are python imports like this? 1st option
    # does not work when running this file by itself. Second doesn't work when importing the file........................................................................................................................



def type_phenotype_step(ctype, this_type_dicts):
    if not this_type_dicts:
        return ''
    empty = ''
    any_pheno = False
    full = f"\t\t\tfor cell in self.cell_list_by_type(self.{ctype.upper()}):\n"
    for cell_dict in this_type_dicts:
        if "phenotypes" in cell_dict.keys():
            any_pheno = True
            full += "\t\t\t\t# WARNING: currently you are responsible for implementing what should happen for each " \
                    "of\n"\
                    "\t\t\t\t# the flags\n"
            full += f"\t\t\t\tchanged_phase, should_be_removed, divides = \\ \n" \
                    f"\t\t\t\t\tcell.dict['current_phenotype'].time_step_phenotype()\n"
            full += f"\t\t\t\tcell.targetVolume = cell.dict['volume_conversion'] * \\ \n" \
                    f"\t\t\t\t\tcell.dict['current_phenotype'].volume.total\n"
    if any_pheno:
        return full
    else:
        return empty


def generate_phenotypes_loops(cell_types, cell_dicts):
    loops = "\n"
    loops += "\t\tif pcp_imp:\n\t\t\tpass\n"
    for ctype in cell_types:
        this_type_dicts = get_dicts_for_type(ctype, cell_dicts)
        loop = type_phenotype_step(ctype, this_type_dicts)
        loops += loop
    return loops


def generate_phenotype_steppable(cell_types, cell_dicts, first=False):
    if first:
        already_imports = False

    else:
        already_imports = True
    loops = generate_phenotypes_loops(cell_types, cell_dicts)

    pheno_step = generate_steppable("Phenotype", 1, True, already_imports=already_imports,
                                    additional_step=loops)
    return pheno_step