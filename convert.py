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

from pathlib import Path

import argparse

from cc3d_xml_gen.gen import make_potts, make_metadata, make_cell_type_plugin, make_cc3d_file, \
    make_contact_plugin, make_diffusion_plug, reconvert_spatial_parameters_with_minimum_cell_volume, make_secretion, \
    reconvert_cell_volume_constraints, decrease_domain, reconvert_time_parameter, make_volume, make_chemotaxis

from cc3d_xml_gen.get_physicell_data import get_cell_constraints, get_secretion_uptake, get_microenvironment, \
    get_dims, get_time, get_chemotaxis
from conversions.secretion import convert_secretion_uptake_data

try:
    from autopep8 import fix_code

    pep_auto = True
except ImportError:
    fix_code = None
    pep_auto = False


def _fix_code(source, options=None, encoding=None, apply_config=False):
    """
    Function to keep the program working if the user doesn't have autopep8. It has the same signature as
    `autopep8.fixcode` but simply returns the 1st argument
    :param source: the string to (not) be converted
    :param options:
    :param encoding:
    :param apply_config:
    :return: source
    """
    return source


if not pep_auto:
    fix_code = _fix_code

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


def default_initial_cell_config(celltypes, xmax, ymax, zmax):
    """
    Returns the default UniformInitializer steppable.

    The default_initial_cell_config function returns an XML string with a configuration for
    the UniformInitializer steppable in CompuCell3D. The configuration is based on the input parameters of celltypes,
    xmax, ymax, and zmax.

    The UniformInitializer steppable is responsible for initializing the initial configuration of cells in the
    simulation. This function sets up a rectangular slab of cells with a specified width and gap size, and restricts
    the cell types to those specified in the celltypes list. The WALL cell type is excluded from the initialization
    process.

    Parameters
    ----------
    celltypes : list of str
        List of cell types.
    xmax : int
        Maximum x dimension of the simulation.
    ymax : int
        Maximum y dimension of the simulation.
    zmax : int
        Maximum z dimension of the simulation.

    Returns
    -------
    str
        Configured XML string.
    """
    beg = '''<Steppable Type="UniformInitializer">
\t<!-- Initial layout of cells in the form of rectangular slab -->
\t<!-- PhysiCell has many complex ways of defining the initial arangement -->
\t<!-- of cells. By default the translator uses a simple configuration, -->
\t<!-- you are responsible for analysing the initialization of the original -->
\t<!-- model and reimplement it accordingly -->
\t<Region>\n'''
    if (zmax != 1 and zmax != 0) and (xmax > 10 and ymax > 10 and zmax > 10):
        box_min = f'\t\t<BoxMin x="{10}" y="{10}" z="{10}"/>\n'
        box_max = f'\t\t<BoxMax x="{xmax - 10}" y="{ymax - 10}" z="{zmax - 10}"/>\n'
    elif zmax != 1 and zmax != 0:
        box_min = f'\t\t<BoxMin x="{1}" y="{1}" z="{1}"/>\n'
        box_max = f'\t\t<BoxMax x="{xmax - 1}" y="{ymax - 1}" z="{zmax - 1}"/>\n'
    else:
        box_min = f'\t\t<BoxMin x="{10}" y="{10}" z="{0}"/>\n'
        box_max = f'\t\t<BoxMax x="{xmax - 10}" y="{ymax - 10}" z="{1}"/>\n'

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
    steppable_string = beg + box_min + box_max + gap + types + end
    return steppable_string


def main(path_to_xml, out_directory=None, minimum_volume=8, max_volume=150 ** 3, name=None):
    """
    Converts a PhysiCell simulation XML into a CompuCell3D simulation folder

    This function is responsible for calling all the methods that will parse PhysiCell's XML and generate the relevant
    CC3D simulation files (steppable file, XML, and main python file). It creates the output folder if one is given and
    places the converted simulation there, if no output path is given the converted simulation will be placed in the
    same directory as the original PhysiCell XML.

    :param max_volume:
    :type max_volume:
    :param minimum_volume:
    :type minimum_volume:
    :param path_to_xml: string for the path to the PhysiCell XML simulation file
    :param out_directory: string path to the output folder
    :param name:
    :return: None
    """

    if minimum_volume is None:
        minimum_volume = 8

    if max_volume is None:
        max_volume = 150 ** 3

    read_before_run = "******************************************************************************************\n" \
                      "PLEASE READ BEFORE RUNNING:\n" \
                      "The translator program does a lot of the conversion process of a Physicell simulation into\n" \
                      "a CompuCell3D simulation. However, just as in CC3D, the user can define several custom " \
                      "modules\n" \
                      "functions, and initial conditions for their PhysiCell model. Translating custom c++ code is\n" \
                      "beyond what an automated program can do (AI-ML language models not withstanding). On top\n" \
                      "of that fact there are several PhysiCell concepts that do not exist in CC3D, and vice-versa.\n" \
                      "Therefore: 1. please read *all comments* the script has placed in the converted files.\n" \
                      "2. You are responsible for creating the initial conditions. Some are defined in csv files\n" \
                      "others are creating with c++.\n" \
                      "3. You are responsible for finding where <user_parameters> is used in the Physicell model,\n" \
                      "and using it in the CC3D model\n" \
                      "4. You are responsible for using chemical field data (e.g., chemotaxis)\n" \
                      "5. You are responsible for defining the use of custom_data for each cell type. Find where \n" \
                      "it is used in PhysiCell and define its use in the CC3D simulation\n" \
                      "******************************************************************************************\n"

    # finding file, making folders
    path_to_xml = Path(path_to_xml)

    xml_dir = path_to_xml.parent

    sim_name = path_to_xml.name.split(".")[0]

    if out_directory is None:
        out_directory = xml_dir.joinpath("CC3D_converted_sim", sim_name)
    else:
        out_directory = Path(out_directory)
    print(f"Creating {out_directory}")
    if not out_directory.exists():
        out_directory.mkdir(parents=True)

    print(f"Creating {out_directory}/Simulation")
    sim_dir = out_directory.joinpath("Simulation")
    if not sim_dir.exists():
        sim_dir.mkdir(parents=True)

    # getting simulation name
    if name is None:
        print("Finding simulation name")
        name = path_to_xml.name.split(".")[0]

    print(f"Creating {out_directory}/{name}.cc3d")

    with open(f"{out_directory}/{name}.cc3d", "w+") as f:
        cc3d, xml_name, main_py_name, steppables_py_name = make_cc3d_file(name=name)
        f.write(cc3d)

    print(f"Loading {path_to_xml}")
    with open(path_to_xml, 'r') as f:
        xml_raw = f.read()

    pcdict = x2d.parse(xml_raw)['PhysiCell_settings']

    print("Generating <Metadata/>")
    metadata_str, n_threads = make_metadata(pcdict)

    print("Extracting space and time data")
    pcdims, ccdims = get_dims(pcdict)
    pctime, cctime = get_time(pcdict)

    print("Generating <Plugin CellType/>")
    ct_str, wall, cell_types, = make_cell_type_plugin(pcdict)

    print("Detecting if the cells are too small or the simulation is too big")
    constraints, any_below, pixel_volumes, minimum_volume = \
        get_cell_constraints(pcdict, ccdims[4], minimum_volume=minimum_volume)
    old_cons = constraints
    old_ccdims = ccdims
    if any_below:
        ccdims, constraints = \
            reconvert_spatial_parameters_with_minimum_cell_volume(constraints, ccdims, pixel_volumes, minimum_volume)
    else:
        constraints = reconvert_cell_volume_constraints(constraints, 1, minimum_volume)

    ccdims, was_above = decrease_domain(ccdims, max_volume=max_volume)

    print("parsing micro environment")
    d_elements = get_microenvironment(pcdict, ccdims[4], pcdims[3], cctime[2], pctime[1])

    d_elements, cctime = reconvert_time_parameter(d_elements, cctime)

    print("Generating <Potts/>")
    potts_str = make_potts(pcdims, ccdims, pctime, cctime)

    with open(os.path.join(sim_dir, "extra_definitions.py"), 'w+') as f:

        f.write(fix_code("cell_constraints=" + str(constraints) + "\n",
                         options={"aggressive": 1}))

    print("Generating <Plugin Contact/>")
    contact_plug = make_contact_plugin(cell_types)

    intializer_step = default_initial_cell_config(cell_types, ccdims[0], ccdims[1], ccdims[2])

    print("Generating diffusion plugin")
    diffusion_string = make_diffusion_plug(d_elements, cell_types, ccdims[6])

    print("Parsing secretion data")
    secretion_uptake_dict = get_secretion_uptake(pcdict)

    conv_sec = convert_secretion_uptake_data(secretion_uptake_dict, cctime[2], pctime[1])

    print("Parsing chemotaxis data")
    taxis_dict = get_chemotaxis(pcdict)

    chemotaxis_plug = make_chemotaxis(taxis_dict)

    print("Generating constraint steppable")
    constraint_step = steppable_gen.generate_constraint_steppable(cell_types, [constraints,
                                                                               conv_sec, taxis_dict], wall,
                                                                  user_data=pcdict["user_parameters"])

    print("Generating secretion steppable")
    secretion_step = steppable_gen.generate_secretion_uptake_step(cell_types, secretion_uptake_dict)

    secretion_plug = make_secretion(secretion_uptake_dict)

    print("Generating phenotype steppable")
    pheno_step = steppable_gen.generate_phenotype_steppable(cell_types, [constraints,
                                                                         conv_sec])

    print("Generating CC3DML")
    cc3dml = "<CompuCell3D>\n"
    cc3dml += "<!--\n" + read_before_run + "-->\n"
    cc3dml += metadata_str + potts_str + ct_str + make_volume() + contact_plug + diffusion_string + secretion_plug + \
              '\n' + chemotaxis_plug + '\n' + \
              intializer_step + "\n\n" + "\n</CompuCell3D>\n"

    print(f"Creating {out_directory}/Simulation/{xml_name}")
    with open(os.path.join(out_directory, f"Simulation/{xml_name}"), "w+") as f:
        f.write(cc3dml)

    print("Merging steppables")

    all_step = '"""\n' + read_before_run + '"""\n' + constraint_step + "\n" + secretion_step + "\n" + pheno_step

    step_names = steppable_gen.get_steppables_names(all_step)

    print("Generating steppables file")

    steppable_gen.generate_steppable_file(sim_dir, f"{steppables_py_name}", fix_code(all_step,
                                                                                     options={"aggressive": 1}))

    print("Generating steppable registration file")
    steppable_gen.generate_main_python(sim_dir, f"{main_py_name}", f"{steppables_py_name}", step_names, read_before_run)

    print("______________\nDONE!!")
    return


parser = argparse.ArgumentParser(description="Converts a Physicell XML file into CompuCell3D .cc3d, .xml, main.py, and"
                                             "steppables.py simulation configuration files.")
parser.add_argument("input", type=str, help="Path to your input PhysiCell XML configuration file")
parser.add_argument("-c", "--cellvolume", type=int, help="(optional) minimum volume the converted cells are allowed to "
                                                         "have (in pixels)", default=None)
parser.add_argument("-v", "--simulationvolume", type=int, help="(optional) maximum volume the CC3D simulation can have",
                    default=None)
parser.add_argument("-o", "--output", help="(optional) output path for the converted files",
                    default=None)
args = parser.parse_args()
main(args.input, out_directory=args.output, minimum_volume=args.cellvolume, max_volume=args.simulationvolume)

