try:
    from .gen_functions import generate_steppable, steppable_imports
except:
    from gen_functions import generate_steppable, steppable_imports  # why are python imports like this? 1st option
    # does not work when running this file by itself. Second doesn't work when importing the file........................................................................................................................


def _apply_volume_constraint(cdict):
    cstr = f'\t\t\tcell.targetVolume = {round(cdict["volume (pixels)"])}'
    cstr += '\n\t\t\t# NOTE: PC does not ' \
            f'have an equivalent parameter, you have to adjust it:\n' \
            f'\t\t\tcell.lambdaVolume = 16\n'
    return cstr


def _apply_surface_constraint(cdict):
    cstr = f'\t\t\tcell.targetSurface = {cdict["surface (pixels)"]}'
    cstr += '\n\t\t\t# NOTE: PC does not ' \
            'have an equivalent parameter, you have to adjust it:'
    cstr += '\n\t\t\tcell.lambdaSurface = 8\n'
    return cstr


# noinspection PyPep8Naming
def apply_CC3D_constraint(cname, cdict):
    if cname.upper() == "VOLUME":
        return _apply_volume_constraint(cdict)
    elif cname.upper() == "SURFACE":
        return _apply_surface_constraint(cdict)
    return ''


# def apply_phenotype()


def cell_type_constraint(ctype, this_type_dicts):
    if not this_type_dicts:
        return ''
    loop = f"\t\tfor cell in self.cell_list_by_type(self.{ctype.upper()}):\n"
    full = loop
    for cell_dict in this_type_dicts:
        if "chemotaxis_dict" in cell_dict.keys():
            line = ""
            for key, value in cell_dict.items():
                if key == "chemotaxis_dict":
                    continue
                field_name = key
                line += f'\t\t\tcd = self.chemotaxisPlugin.addChemotaxisData(cell, "{field_name}")\n'
                line += f'\t\t\tcd.setLambda({value}*100)\n'

            full += line
            continue
        for key, value in cell_dict.items():
            if key == "phenotypes" and bool(cell_dict["phenotypes"]):
                line = "\t\t\tif pcp_imp:\n"
                line += f"\t\t\t\tcell.dict['{key}']=self.phenotypes['{ctype}']\n"
                line += f"\t\t\t\tcell.dict['current_phenotype'] = cell.dict['{key}']" \
                        f"['{cell_dict['phenotypes_names'][0]}'].copy()\n"
                line += f"\t\t\t\tcell.dict['volume_conversion'] = cell.targetVolume / \\\n" \
                        f"\t\t\t\t\tcell.dict['current_phenotype'].current_phase.volume.total\n"
            elif key == "custom_data":
                line = f"\t\t\t# NOTE: you are responsible for finding how this data" \
                       f"is used in the original model\n\t\t\t# and re-implementing in CC3D" \
                       f"\n\t\t\tcell.dict['{key}']={value}\n"
            elif type(value) == str:
                line = f"\t\t\tcell.dict['{key}']='{value}'\n"

            elif type(value) == list:
                line = f"\t\t\tcell.dict['{key}']={value}\n"
            elif type(value) == dict and key != "phenotypes":
                clean_value = value.copy()
                to_pop = []
                for subkey in value.keys():

                    if "comment" in subkey:
                        to_pop.append(subkey)
                [clean_value.pop(p) for p in to_pop]
                line = f"\t\t\tcell.dict['{key}']={clean_value}\n"
            # else:
            #     line = f"\t\t\tcell.dict['{key}']={value}\n"
            if key in ["volume", "surface"]:
                line += apply_CC3D_constraint(key, value)
            full += line
    return full + '\n\n'


def get_dicts_for_type(ctype, cell_dicts):
    ds = []
    for d in cell_dicts:
        if ctype in d.keys():
            ds.append(d[ctype])
    return ds


def generate_constraint_loops(cell_types, cell_dicts):
    loops = "\n"
    for ctype in cell_types:
        this_type_dicts = get_dicts_for_type(ctype, cell_dicts)
        loop = cell_type_constraint(ctype, this_type_dicts)
        loops += loop
    return loops


def initialize_phenotypes(constraint_dict):
    pheno_str = "\n\t\tif pcp_imp:\n"
    pheno_str += "\t\t\tself.phenotypes = {}\n"
    for ctype, cdict in constraint_dict.items():
        if "phenotypes" in cdict.keys():

            pheno_str += f"\t\t\tdt = 1/self.mcs_to_time\n"
            pheno_str += f"\t\t\tself.phenotypes['{ctype}']" + "= {}\n"
            for phenotype, _data in cdict["phenotypes"].items():
                data = _data
                if data is None:
                    data = {}
                time_unit = "None"
                if _data is not None and "rate units" in _data.keys():
                    time_unit = _data["rate units"].split("/")[-1]
                fixed = []
                duration = []
                if _data is not None:
                    for fix, dur in _data["phase durations"]:
                        duration.append(dur)
                        if fix == "TRUE":
                            ff = True
                        else:
                            ff = False
                        fixed.append(ff)
                else:
                    fixed.append(False)
                    duration.append(None)
                nuclear_fluid = []
                nuclear_solid = []
                cyto_fluid = []
                cyto_solid = []
                cyto_to_nucl = []
                if _data is None:
                    data['fluid fraction'] = [.75]
                elif 'fluid fraction' not in _data.keys():
                    data['fluid fraction'] = [.75] * len(_data["phase durations"])
                # if 'fluid fraction' not in data.keys():
                #     data['fluid fraction'] = [.75]*len(data["phase durations"])
                if _data is not None and 'fluid fraction' in _data.keys() and 'nuclear volume' in _data.keys() and \
                        'total' in _data.keys():
                    for fluid, nucl, total in zip(_data['fluid fraction'], _data['nuclear volume'], _data['total']):
                        nfl = fluid * nucl
                        nuclear_fluid.append(nfl)
                        nuclear_solid.append(nucl - nfl)
                        cytt = total - nucl
                        cytf = fluid * cytt
                        cyts = cytt - cytf
                        cyto_fluid.append(cytf)
                        cyto_solid.append(cyts)
                        cyto_to_nucl.append(cytt / (1e-16 + nucl))
                else:
                    if _data is not None:
                        nuclear_fluid = [None] * len(_data["phase durations"])
                        nuclear_solid = [None] * len(_data["phase durations"])
                        cyto_fluid = [None] * len(_data["phase durations"])
                        cyto_solid = [None] * len(_data["phase durations"])
                        cyto_to_nucl = [None] * len(_data["phase durations"])
                    else:
                        nuclear_fluid = [None]
                        nuclear_solid = [None]
                        cyto_fluid = [None]
                        cyto_solid = [None]
                        cyto_to_nucl = [None]
                if _data is None:
                    data['calcified fraction'] = [0]
                elif 'calcified fraction' not in _data.keys():
                    data['calcified fraction'] = [0] * len(_data["phase durations"])

                if _data is None:
                    data['cytoplasm biomass change rate'] = [None]
                elif 'cytoplasm biomass change rate' not in _data.keys():
                    data['cytoplasm biomass change rate'] = [None] * len(_data["phase durations"])

                if _data is None:
                    data['nuclear biomass change rate'] = [None]
                elif 'nuclear biomass change rate' not in _data.keys():
                    data['nuclear biomass change rate'] = [None] * len(_data["phase durations"])
                if _data is None:
                    data['calcification rate'] = [None]
                elif 'calcification rate' not in _data.keys():
                    data['calcification rate'] = [None] * len(_data["phase durations"])

                if _data is None:
                    data['fluid change rate'] = [None]
                elif 'fluid change rate' not in _data.keys():
                    data['fluid change rate'] = [None] * len(_data["phase durations"])

                pheno_str += f"\t\t\tphenotype = pcp.get_phenotype_by_name('{phenotype}')\n"
                pheno_str += f"\t\t\tself.phenotypes['{ctype}']['{phenotype}'] = phenotype(dt=dt, \n\t\t\t\t" \
                             f"time_unit='{time_unit}', \n\t\t\t\tfixed_durations={fixed},  " \
                             f"\n\t\t\t\tphase_durations={duration}, \n\t\t\t\t" \
                             f"cytoplasm_volume_change_rate={data['cytoplasm biomass change rate']}, \n\t\t\t\t" \
                             f"nuclear_volume_change_rate={data['nuclear biomass change rate']}, \n\t\t\t\t" \
                             f"calcification_rate={data['calcification rate']}, \n\t\t\t\t" \
                             f"calcified_fraction={data['calcified fraction']}, \n\t\t\t\t" \
                             f"target_fluid_fraction={data['fluid fraction']}, \n\t\t\t\t" \
                             f"nuclear_fluid={nuclear_fluid}, \n\t\t\t\t" \
                             f"nuclear_solid={nuclear_solid}, \n\t\t\t\t" \
                             f"nuclear_solid_target={nuclear_solid}, \n\t\t\t\t" \
                             f"cytoplasm_fluid={cyto_fluid}, \n\t\t\t\t" \
                             f"cytoplasm_solid={cyto_solid}, \n\t\t\t\t" \
                             f"cytoplasm_solid_target={cyto_solid}, \n\t\t\t\t" \
                             f"target_cytoplasm_to_nuclear_ratio={cyto_to_nucl}, \n\t\t\t\t" \
                             f"fluid_change_rate={data['fluid change rate']})\n"

    return pheno_str


def generate_constraint_steppable(cell_types, cell_type_dicts, wall, first=True, user_data=""):
    already_imports = not first
    loops = generate_constraint_loops(cell_types, cell_type_dicts)
    if not wall:
        wall_str = "\t\tself.shared_steppable_vars['constraints'] = self"
    else:
        wall_str = "\t\tself.build_wall(self.WALL)\n\t\tself.shared_steppable_vars['constraints'] = self"
    pheno_init = initialize_phenotypes(cell_type_dicts[0])
    constraint_step = generate_steppable("Constraints", 1, False, minimal=True, already_imports=already_imports,
                                         additional_start=pheno_init + loops + wall_str, user_data=user_data)
    return constraint_step


if __name__ == "__main__":
    dicts = [{"a": {"exemplo": 1, "asdasd": 3}, "b": {"kkkkkk": "popop"}},
             {"a": {"exeo": 1, "uip": 3}, "b": {"oooooo": "lololo"}}]
    print(generate_constraint_steppable(["a", "b"], dicts, False))
