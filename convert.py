# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 09:45:52 2022

@author: Juliano Ferrari Gianlupi
"""

# import string
# import copy
import os
import xmltodict as x2d
import steppable_gen

from cc3d_xml_gen.gen import make_potts, make_metadata, make_cell_type_plugin, make_cc3d_file, \
    make_contact_plugin, make_diffusion_plug, reconvert_spatial_parameters_with_minimum_cell_volume, make_secretion

from cc3d_xml_gen.get_physicell_data import get_cell_constraints, get_secretion, get_microenvironment, get_dims, \
    get_time
from conversions.secretion import convert_secretion_data

# TODO:
#  - phenotypes
#  - bash script
#  - cell configuration
#  - mechanics (e.g. elasticity)
#  - dynamic simulation name gen; custom sim name gen
#  - cc3dml file has to always have the correct file names


# defines conversion factors to meter
_space_convs = {"micron": 1e-6,
                "micrometer": 1e-6,
                "micro": 1e-6,
                "milli": 1e-3,
                "millimeter": 1e-3,
                "nano": 1e-9,
                "nanometer": 1e-9,
                'meter': 1
                }

# defines conversion factors to minutes
_time_convs = {"millisecond": 1e-3 / 60,
               "milliseconds": 1e-3 / 60,
               "microsecond": 1e-6 / 60,
               "microseconds": 1e-6 / 60,
               "second": 1 / 60,
               "s": 1 / 60,
               "seconds": 1 / 60,
               "hours": 60,
               "hour": 60,
               "h": 60,
               "day": 24 * 60,
               "days": 24 * 60,
               "week": 7 * 24 * 60,
               "weeks": 7 * 24 * 60,
               "minutes": 1,
               "min": 1}


def extra_for_testing(celltypes, xmax, ymax, zmax):
    beg = '''<Steppable Type="UniformInitializer">
\t<!-- Initial layout of cells in the form of rectangular slab -->
\t<Region>\n'''

    box_min = f'\t\t<BoxMin x="{max(1, xmax - 50)}" y="{max(1, ymax - 50)}" z="{max(1, xmax - 50)}"/>\n'
    box_max = f'\t\t<BoxMax x="{xmax - 10}" y="{ymax - 10}" z="{zmax - 10}"/>\n'

    gap = "\t\t<Gap>0</Gap>\n\t\t<Width>7</Width>\n"

    types = ''
    for t in celltypes:
        if t.upper() == "WALL":
            continue
        types += f"{t},"

    types = types[:-1]

    types = "\t\t<Types>" + types + "</Types>\n"

    end = '''\t</Region>
</Steppable>'''

    return beg + box_min + box_max + gap + types + end


if __name__ == "__main__":

    print("Running test")

    print("Creating ./test")
    out_folder = r"./test"
    if not os.path.isdir(out_folder):
        os.mkdir(out_folder)

    print("Creating ./test/Simulation")

    out_sim_f = os.path.join(out_folder, "Simulation")
    if not os.path.isdir(out_sim_f):
        os.mkdir(out_sim_f)

    print("Creating test/test.cc3d")

    with open("test/test.cc3d", "w+") as f:
        f.write(make_cc3d_file())

    # todo: remember to extract the name of the original PC sim and use it for naming the cc3d sim
    example_path = r"./example_pcxml/" + "annotated_cancer_immune3D_flat.xml"
    # example_path = r"./example_pcxml/" + "virus_macrophage_flat.xml"

    print(f"Loading {example_path}")
    with open(example_path, 'r') as f:
        xml_raw = f.read()
    pcdict = x2d.parse(xml_raw)['PhysiCell_settings']

    print("Generating <Metadata/>")
    metadata_str, n_threads = make_metadata(pcdict)

    print("Extracting space and time data")
    pcdims, ccdims = get_dims(pcdict)
    pctime, cctime = get_time(pcdict)


    print("Generating <Plugin CellType/>")
    ct_str, wall, cell_types, = make_cell_type_plugin(pcdict)

    constraints, any_below, pixel_volumes, minimum_volume = \
        get_cell_constraints(pcdict, ccdims[4], minimum_volume=8)
    old_cons = constraints
    old_ccdims = ccdims
    if any_below:
        ccdims, constraints = \
            reconvert_spatial_parameters_with_minimum_cell_volume(constraints, ccdims, pixel_volumes, minimum_volume)
    print("Generating <Potts/>")
    potts_str = make_potts(pcdims, ccdims, pctime, cctime)

    with open(os.path.join(out_sim_f, "extra_definitions.py"), 'w+') as f:
        f.write("cell_constraints=" + str(constraints) + "\n")

    print("Generating <Plugin Contact/>")
    contact_plug = make_contact_plugin(cell_types)

    extra = extra_for_testing(cell_types, ccdims[0], ccdims[1], ccdims[2])

    print("parsing micro environment")
    d_elements = get_microenvironment(pcdict, ccdims[4], pcdims[3], cctime[2], pctime[1])

    print("Generating diffusion plugin")
    diffusion_string = make_diffusion_plug(d_elements, cell_types, False)

    print("Parsing secretion data")
    secretion_dict = get_secretion(pcdict)

    conv_sec = convert_secretion_data(secretion_dict, cctime[2], pctime[1])

    print("Generating constraint steppable")
    constraint_step = steppable_gen.generate_constraint_steppable(cell_types, [constraints,
                                                                               conv_sec], wall)

    print("Generating secretion steppable")
    secretion_step = steppable_gen.generate_secretion_step(cell_types, secretion_dict)

    secretion_plug = make_secretion(secretion_dict)

    print("Generating CC3DML")
    cc3dml = "<CompuCell3D>\n"
    cc3dml += metadata_str + potts_str + ct_str + contact_plug + diffusion_string + secretion_plug + '\n' + extra + \
              "\n</CompuCell3D>\n"

    print(f"Creating {out_sim_f}/test.xml")
    with open(os.path.join(out_sim_f, "test.xml"), "w+") as f:
        f.write(cc3dml)

    print("Merging steppables") 

    all_step = constraint_step + "\n" + secretion_step

    step_names = steppable_gen.get_steppables_names(all_step)

    print("Generating steppables file")

    steppable_gen.generate_steppable_file(out_sim_f, "steppable_test.py", all_step)

    print("Generating steppable registration file")
    steppable_gen.generate_main_python(out_sim_f, "main_test.py", "steppable_test.py", step_names)

    print("______________\nDONE!!")
