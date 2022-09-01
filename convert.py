# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 09:45:52 2022

@author: Juliano Ferrari Gianlupi
"""

import sys
# import string
# import copy
import os
import warnings
import shutil as sh
import xmltodict as x2d
import steppable_gen

from itertools import combinations

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


def get_dims(pcdict, space_convs=_space_convs):
    xmin = float(pcdict['domain']['x_min']) if "x_min" in pcdict['domain'].keys() else None
    xmax = float(pcdict['domain']['x_max']) if "x_max" in pcdict['domain'].keys() else None

    ymin = float(pcdict['domain']['y_min']) if "y_min" in pcdict['domain'].keys() else None
    ymax = float(pcdict['domain']['y_max']) if "y_max" in pcdict['domain'].keys() else None

    zmin = float(pcdict['domain']['z_min']) if "z_min" in pcdict['domain'].keys() else None
    zmax = float(pcdict['domain']['z_max']) if "z_max" in pcdict['domain'].keys() else None

    units = pcdict['overall']['space_units'] if 'overall' in pcdict.keys() and \
                                                'space_units' in pcdict['overall'].keys() else 'micron'

    autoconvert_space = True
    if units not in space_convs.keys():
        message = f"WARNING: {units} is not part of known space units. Automatic space-unit conversion disabled." \
                  f"Available units for auto-conversion are:\n{space_convs.keys()}"
        warnings.warn(message)
        autoconvert_space = False

    # the dx/dy/dz tags mean that for every voxel there are dx space-units.
    # therefore [dx] = [space-unit/voxel]. Source: John Metzcar

    dx = float(pcdict['domain']['dx']) if "dx" in pcdict['domain'].keys() else 1
    dy = float(pcdict['domain']['dy']) if "dy" in pcdict['domain'].keys() else 1
    dz = float(pcdict['domain']['dz']) if "dz" in pcdict['domain'].keys() else 1

    # print(dx, dy, dz, type(dx), type(dy), type(dz), )
    if not dx == dy == dz:
        message = "WARNING! Physicell's dx/dy/dz are not all the same: " \
                  f"dx={dx}, dy={dy}, dz={dz}\n" \
                  f"Using {max(min(dx, dy, dz), 1)}"
        warnings.warn(message)

    ds = max(min(dx, dy, dz), 1)

    diffx = 1 if xmin is None or xmax is None else round(xmax - xmin)
    diffy = 1 if ymin is None or ymax is None else round(ymax - ymin)
    diffz = 1 if zmin is None or zmax is None else round(ymax - zmin)

    cc3dx = round(diffx / ds)
    cc3dy = round(diffy / ds)
    cc3dz = round(diffz / ds)

    cc3dds = cc3dx / diffx  # pixel/unit

    # [cc3dds] = pixel/unit
    # [cc3dds] * unit = pixel

    cc3dspaceunitstr = f"1 pixel = {cc3dds} {units}"

    return ((xmin, xmax), (ymin, ymax), (zmin, zmax), units), \
           (cc3dx, cc3dy, cc3dz, cc3dspaceunitstr, cc3dds, autoconvert_space)


def get_time(pcdict, time_convs=_time_convs):
    mt = float(pcdict['overall']['max_time']['#text']) if "max_time" in pcdict['overall'].keys() and \
                                                          '#text' in pcdict['overall']['max_time'].keys() else 100000

    mtunit = pcdict['overall']['max_time']['@units'] if "max_time" in pcdict['overall'].keys() and '@units' in \
                                                        pcdict['overall']['max_time'].keys() else None

    time_unit = pcdict['overall']['time_units'] if "time_units" in pcdict['overall'].keys() else None

    if mtunit != time_unit:
        message = f"Warning: Psysicell time units in " \
                  "\n`<overall>\n\t<max_time units=...`\ndiffers from\n" \
                  f"`<time_units>unit</time_units>`.\nUsing: {time_unit}"
        warnings.warn(message)

    autoconvert_time = True
    if time_unit not in time_convs.keys():
        message = f"WARNING: {time_unit} is not part of known time units. Automatic time-unit conversion disabled." \
                  f"Available units for auto-conversion are:\n{time_convs.keys()}"
        warnings.warn(message)
        autoconvert_time = False

    mechdt = float(pcdict['overall']['dt_mechanics']['#text']) if "dt_mechanics" in pcdict['overall'].keys() else None

    steps = round(mt / mechdt)

    cc3ddt = 1 / (mt / steps)  # MCS/unit

    # [cc3ddt] = MCS/unit
    # [cc3ddt] * unit = MCS

    cc3dtimeunitstr = f"1 MCS = {cc3ddt} {time_unit}"

    # timeconvfact = 1/cc3ddt

    return (mt, time_unit, mechdt), (steps, cc3dtimeunitstr, cc3ddt, autoconvert_time)


def get_parallel(pcdict):
    return int(pcdict['parallel']['omp_num_threads']) if 'parallel' in pcdict.keys() and \
                                                         'omp_num_threads' in pcdict['parallel'].keys() else 1


def make_potts(pcdict):
    # todo: figure out the spatial dimensions. What dx/dy/dz mean in
    #     <domain>
    # 		<x_min>-400</x_min>
    # 		<x_max>400</x_max>
    # 		<y_min>-400</y_min>
    # 		<y_max>400</y_max>
    # 		<z_min>-10</z_min>
    # 		<z_max>10</z_max>
    # 		<dx>20</dx>
    # 		<dy>20</dy>
    # 		<dz>20</dz>
    # 		<use_2D>true</use_2D>
    # 	</domain>
    pcdims, ccdims = get_dims(pcdict)

    # space_units_str = f'"1 pixel = 1 {space_units}"'

    pctime, cctime = get_time(pcdict)

    # still need to implement space units
    potts_str = f""" 
<Potts>
   <!-- Basic properties of CPM (GGH) algorithm -->
   <Space_Units>{ccdims[3]}</Space_Units>
   <Pixel_to_Space units="pixel/{pcdims[3]}">{ccdims[4]}</Pixel_to_Space>
   <Dimensions x="{ccdims[0]}" y="{ccdims[1]}" z="{ccdims[2]}"/>
   <Time_Units>{cctime[1]}</Time_Units>
   <MCS_to_Time units="MCS/{pctime[1]}">{cctime[2]}</MCS_to_Time>
   <Steps>{cctime[0]}</Steps>
   <!-- As the frameworks of CC3D and PhysiCell are very different -->
   <!-- PC doesn't have some concepts that CC3D does. Temperature is one of -->
   <!-- them, so the translation script leaves its tunning as an exercise-->
   <!-- for the reader -->
   <Temperature>10.0</Temperature>
   <!-- Same deal for neighbor order as for temperature-->
   <NeighborOrder>1</NeighborOrder>
   <!-- <Boundary_x>Periodic</Boundary_x> -->
   <!-- <Boundary_y>Periodic</Boundary_y> -->
</Potts>\n"""

    return potts_str, pcdims, ccdims, pctime, cctime


def make_metadata(pcdict, out=100):
    threads = get_parallel(pcdict)

    metadata = f'''
<Metadata>
  <!-- Basic properties simulation -->
  <NumberOfProcessors>{threads}</NumberOfProcessors>
  <DebugOutputFrequency>{out}</DebugOutputFrequency>
  <!-- <NonParallelModule Name="Potts"/> -->
</Metadata>\n'''

    return metadata, threads


def get_boundary_wall(pcdict):
    if "options" in pcdict.keys() and "virtual_wall_at_domain_edge" in pcdict['options'].keys():
        wall = pcdict["options"]["virtual_wall_at_domain_edge"].upper()
        if wall == "TRUE":
            return True
        else:
            return False
    else:
        return False
    return False


def get_cell_volume(subdict):
    if 'phenotype' in subdict.keys() and 'volume' in subdict['phenotype'].keys() and \
            'total' in subdict['phenotype']['volume'].keys():
        volume = float(subdict['phenotype']['volume']['total']['#text'])
        units = subdict['phenotype']['volume']['total']['@units']
        return volume, units
    return None, None


def make_cell_type_tags(pcdict):
    s = ''
    cell_types = []
    idx = 1

    # volumes = {}

    for child in pcdict['cell_definitions']['cell_definition']:
        # print(child.tag, child.attrib, child.text)
        name = child['@name'].replace(" ", "_")
        cell_types.append(name)
        ctt = f'\t<CellType TypeId="{idx}" TypeName="{name}"/>\n'
        s += ctt

        idx += 1

    create_wall = get_boundary_wall(pcdict)

    if create_wall:
        s += f'\t<CellType Freeze="" TypeId="{idx}" TypeName="WALL"/>\n'
        cell_types.append("WALL")

    return s, create_wall, cell_types


def make_cell_type_plugin(pcdict):
    ct_str = '\n<Plugin Name="CellType">\n\t' \
             '<CellType TypeId="0" TypeName="Medium"/>\n'
    typesstr, wall, cell_types = make_cell_type_tags(pcdict)

    ct_str += typesstr
    ct_str += '</Plugin>'

    return ct_str, wall, cell_types


def get_cell_mechanics(subdict):
    if 'phenotype' in subdict.keys() and 'mechanics' in subdict['phenotype'].keys():
        d = {}
        for key, item in subdict['phenotype']['mechanics'].items():
            if key != "options":
                d[key] = {'units': item['@units'],
                          'value': float(item['#text'])}
    else:
        return None

    return d


def get_cell_constraints(pcdict, space_unit, time_unit):
    constraints = {}

    for child in pcdict['cell_definitions']['cell_definition']:
        ctype = child['@name'].replace(" ", "_")
        constraints[ctype] = {}
        volume, unit = get_cell_volume(child)
        dim = int(unit.split("^")[-1])
        volumepx = volume * (space_unit ** dim)
        constraints[ctype]["volume"] = {f"volume ({unit})": volume,
                                        "volume (pixels)": volumepx}
        constraints[ctype]["mechanics"] = get_cell_mechanics(child)

    return constraints


def make_cc3d_file(name=None):
    if name is None:
        cc3d = '''
<Simulation version="4.3.0">
   <XMLScript Type="XMLScript">Simulation/test.xml</XMLScript>
    <PythonScript Type="PythonScript">Simulation/test.py</PythonScript> 
    <Resource Type="Python">Simulation/testSteppables.py</Resource> 
    <Resource Type="Python">Simulation/extra_definitions.py</Resource> 
</Simulation>\n'''
        return cc3d
    else:
        cc3d = '''
<Simulation version="4.3.0">
   <XMLScript Type="XMLScript">Simulation/{name}.xml</XMLScript>
   <PythonScript Type="PythonScript">Simulation/{name}.py</PythonScript>
   <Resource Type="Python">Simulation/{name}Steppables.py</Resource>
   <Resource Type="Python">Simulation/extra_definitions.py</Resource> 
</Simulation>\n'''
        return cc3d


def make_contact_plugin(celltypes):
    # replace ASAP with contact flex

    combs = list(combinations(celltypes, 2))

    for t in celltypes:
        combs.append((t, t))

    combs.reverse()

    contact_plug = """
<Plugin Name="Contact">
\t<!-- PhysiCell doesn't have an equivalent to this plugin. Its  -->
\t<!-- tunning and deciding on the neighbor order is left as an -->
\t<!-- exerise to the reader. -->
\t<!-- A better option (to be implemented) is to use the adhesion flex -->
\t<!-- Specification of adhesion energies -->
\t<Energy Type1="Medium" Type2="Medium">10.0</Energy>\n"""

    # 1 make the medium contact energies
    me = ""
    for t in celltypes:
        me += f'\t<Energy Type1="Medium" Type2="{t}">10.0</Energy>\n'

    # 2 make the combination energies

    ce = ""
    for t1, t2 in combs:
        ce += f'\t<Energy Type1="{t1}" Type2="{t2}">10.0</Energy>\n'

    contact_plug += me + ce + "\t<NeighborOrder>3</NeighborOrder>\n</Plugin>"
    return contact_plug


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


def get_space_time_from_diffusion(unit):
    parts = unit.split("/")
    timeunit = parts[-1]
    spaceunit = parts[0].split("^")[0]
    return spaceunit, timeunit


def get_secretion(pcdict):
    # will have to be done in python
    sec_data = {}
    for child in pcdict['cell_definitions']['cell_definition']:
        ctype =child['@name'].replace(" ", "_")
        sec_data[ctype] = {}
        sec_list = child['phenotype']['secretion']['substrate']
        for sec in sec_list:
            substrate = sec["@name"].replace(" ", "_")
            sec_data[ctype][substrate] = {}
            sec_data[ctype][substrate]['secretion_rate'] = float(sec['secretion_rate']['#text'])
            sec_data[ctype][substrate]['secretion_unit'] = sec['secretion_rate']['@units']
            sec_data[ctype][substrate]['secretion_target'] = float(sec['secretion_target']['#text'])
            sec_data[ctype][substrate]['uptake_rate'] = float(sec['uptake_rate']['#text'])
            sec_data[ctype][substrate]['uptake_unit'] = sec['uptake_rate']['@units']
            sec_data[ctype][substrate]['net_export'] = float(sec['net_export_rate']['#text'])
            sec_data[ctype][substrate]['net_export_unit'] = sec['net_export_rate']['@units']
    return sec_data



def get_microenvironment(pcdict, space_factor, space_unit, time_factor, time_unit, autoconvert_time=True,
                         autoconvert_space=True, space_convs=_space_convs, time_convs=_time_convs):
    diffusing_elements = {}
    fields = pcdict['microenvironment_setup']['variable']
    for subel in fields:

        auto_s_this = autoconvert_space
        auto_t_this = autoconvert_time

        diffusing_elements[subel['@name']] = {}
        diffusing_elements[subel['@name']]["concentration_units"] = subel["@units"]
        diffusing_elements[subel['@name']]["D_w_units"] = \
            float(subel['physical_parameter_set']['diffusion_coefficient']['#text'])

        this_space, this_time = get_space_time_from_diffusion(subel['physical_parameter_set']['diffusion_coefficient']
                                                              ['@units'])

        if this_space != space_unit:
            message = f"WARNING: space unit found in diffusion coefficient of {subel['@name']} does not match" \
                      f"space unit found while converting <overall>:\n\t<overall>:{space_unit};\n\t{subel['@name']}:{this_space}" \
                      f"\nautomatic space-unit conversion for {subel['@name']} disabled"
            warnings.warn(message)
            auto_s_this = False
            # space_conv_factor = 1
        if this_time != time_unit:
            message = f"WARNING: time unit found in diffusion coefficient of {subel['@name']} does not match" \
                      f"space unit found while converting <overall>:\n\t<overall>:{time_unit};" \
                      f"\n\t{subel['@name']}:{this_time}" \
                      f"\nautomatic time-unit conversion for {subel['@name']} disabled"
            warnings.warn(message)
            auto_t_this = False
            # time_conv_factor = 1

        diffusing_elements[subel['@name']]["auto"] = (auto_s_this, auto_t_this)

        if auto_s_this:
            space_conv_factor = space_factor
        else:
            space_conv_factor = 1

        if auto_t_this:
            time_conv_factor = time_factor
        else:
            time_conv_factor = 1

        D = diffusing_elements[subel['@name']]["D_w_units"] * space_conv_factor * space_conv_factor / time_conv_factor
        diffusing_elements[subel['@name']]["D"] = D

        # [cc3dds] = pixel/unit
        # [cc3dds] * unit = pixel -> pixel^2 = ([cc3dds] * unit)^2
        # [cc3ddt] = MCS/unit
        # [cc3ddt] * unit = MCS -> 1/MCS = 1/([cc3ddt] * unit)

        if auto_s_this and auto_t_this:
            diffusing_elements[subel['@name']]["D_conv_factor_text"] = f"1 pixel^2/MCS" \
                                                                       f" = {space_conv_factor * space_conv_factor / time_conv_factor}" \
                                                                       f"{space_unit}^2/{time_unit}"
            diffusing_elements[subel['@name']]["D_conv_factor"] = space_conv_factor * space_conv_factor / \
                                                                  time_conv_factor
            diffusing_elements[subel['@name']]["D_og_unit"] = f"{space_unit}^2/{time_unit}"
        else:
            diffusing_elements[subel['@name']]["D_conv_factor_text"] = "disabled autoconversion"
            diffusing_elements[subel['@name']]["D_conv_factor"] = 1
            diffusing_elements[subel['@name']]["D_og_unit"] = "disabled autoconversion, not known"
        diffusing_elements[subel['@name']]["gamma_w_units"] = float(subel['physical_parameter_set']['decay_rate']
                                                                    ['#text'])
        gamma = diffusing_elements[subel['@name']]["gamma_w_units"] / time_conv_factor
        diffusing_elements[subel['@name']]["gamma"] = gamma
        if auto_t_this:
            diffusing_elements[subel['@name']][
                "gamma_conv_factor_text"] = f"1/MCS = {1 / time_conv_factor} 1/{time_unit}"
            diffusing_elements[subel['@name']]["gamma_conv_factor"] = 1 / time_conv_factor
            diffusing_elements[subel['@name']]["gamma_og_unit"] = f"1/{time_unit}"
        else:
            diffusing_elements[subel['@name']]["gamma_conv_factor_text"] = "disabled autoconversion"
            diffusing_elements[subel['@name']]["gamma_conv_factor"] = 1
            diffusing_elements[subel['@name']]["gamma_og_unit"] = "disabled autoconversion, not known"

        diffusing_elements[subel['@name']]["initial_condition"] = subel['initial_condition']['#text']

        diffusing_elements[subel['@name']]["dirichlet"] = subel['Dirichlet_boundary_condition']['@enabled']
        diffusing_elements[subel['@name']]["dirichlet_value"] = float(subel['Dirichlet_boundary_condition']['#text'])

    return diffusing_elements


def make_diffusion_plug(diffusing_elements, celltypes, flag_2d):
    header = f'\n\n\t<Steppable Type="DiffusionSolverFE">\t<!-- The conversion uses DiffusionSolverFE by default. You may ' \
             f'wish to use another diffusion solver-->'

    full_str = header

    for key, item in diffusing_elements.items():
        name = key.replace(" ", "_")

        # diffusion data
        df_str = f'\t\t<DiffusionField Name="{name}">\n\t\t\t<DiffusionData>\n\t\t\t\t<FieldName>{name}</FieldName>\n'
        conc_units = f'\t\t\t\t<Concentration_units>{item["concentration_units"]}</Concentration_units>\n'
        og_D = f'\t\t\t\t<Original_diffusion_constant D="{item["D_w_units"]}" units= "{item["D_og_unit"]}"/>\n'
        # conv = f'\t\t\t\t<CC3D_to_original units="(pixel^2/MCS)/(item["D_og_unit"])">{item["D_conv_factor"]}' \
        #        '</CC3D_to_original>'
        D_str = f'\t\t\t\t<GlobalDiffusionConstant>{item["D"]}</GlobalDiffusionConstant>\n'
        og_g = f'\t\t\t\t<Original_decay_constant gamma="{item["gamma_w_units"]}" units= "{item["gamma_og_unit"]}"/>\n'
        g_str = f'\t\t\t\t<GlobalDecayConstant>{item["gamma"]}</GlobalDecayConstant>\n'

        init_cond_warn = '\t\t\t\t<!-- CC3D allows for diffusing fields initial conditions, if one was detected it ' \
                         'will -->\n' \
                         '\t\t\t\t<!-- be used here. For several reasons it may not work, if something looks wrong with' \
                         ' -->\n' \
                         '\t\t\t\t<!-- your diffusing field at the start of the simulation this may be the reason.' \
                         ' -->\n' \
                         '\t\t\t\t<!-- CC3D also allows the diffusing field initial condition to be set by a file. ' \
                         'Conversion of a -->\n' \
                         '\t\t\t\t<!-- PhysiCell diffusing field initial condition file into a CC3D compliant one is ' \
                         'left as -->\n' \
                         '\t\t\t\t<!-- an exercise to the reader. -->\n'

        init_cond = f'\t\t\t\t <InitialConcentrationExpression>{item["initial_condition"]}<' \
                    f'/InitialConcentrationExpression>' \
                    '\n\t\t\t\t<!-- <ConcentrationFileName>INITIAL CONCENTRATION FIELD - typically a file with ' \
                    'path Simulation/NAME_OF_THE_FILE.txt</ConcentrationFileName> -->'

        het_warning = "\n\t\t\t\t<!-- CC3D allows the definition of D and gamma on a cell type basis: -->\n"
        cells_str = ""
        for t in celltypes:
            cells_str += f'\t\t\t\t<!--<DiffusionCoefficient CellType="{t}">{item["D"]}</DiffusionCoefficient>-->\n'
            cells_str += f'\t\t\t\t<!--<DecayCoefficient CellType="{t}">{item["gamma"]}</DecayCoefficient>-->\n'
        close_diff_data = "\t\t\t</DiffusionData>\n"

        # boundary conditions

        bc_head = '\t\t\t<BoundaryConditions>\n\t\t\t\t<!-- PhysiCell has either Dirichlet boundary conditions (i.e. ' \
                  'constant ' \
                  'value) -->\n\t\t\t\t<!-- or "free floating" boundary conditions (i.e., constant flux = 0). -->' \
                  '\n\t\t\t\t<!-- CC3D ' \
                  'allows ' \
                  'for more control of boundary conditions, you may want to revisit the issue. -->\n'
        if item['dirichlet']:
            bc_body = f'\t\t\t\t<Plane Axis="X">\n\t\t\t\t\t<ConstantValue PlanePosition="Min" Value=' \
                      f'"{item["dirichlet_value"]}"/>\n\t\t\t\t\t<ConstantValue PlanePosition="Max" Value=' \
                      f'"{item["dirichlet_value"]}"/>\n\t\t\t\t\t<!-- Other options are (examples): -->\n\t\t\t\t\t' \
                      f'<!--<ConstantDerivative PlanePosition="Min" Value="10.0"/> -->\n\t\t\t\t\t<!--' \
                      f'<ConstantDerivative PlanePosition="Max" Value="10.0"/> -->\n\t\t\t\t\t<!--<Periodic/>-->' \
                      '\t\t\t\t</Plane>\n' \
                      f'\t\t\t\t<Plane Axis="Y">\n\t\t\t\t\t<ConstantValue PlanePosition="Min" Value=' \
                      f'"{item["dirichlet_value"]}"/>\n\t\t\t\t\t<ConstantValue PlanePosition="Max" Value=' \
                      f'"{item["dirichlet_value"]}"/>\n\t\t\t\t\t<!-- Other options are (examples): -->\n\t\t\t\t\t' \
                      f'<!--<ConstantDerivative PlanePosition="Min" Value="10.0"/> -->\n\t\t\t\t\t<!--' \
                      f'<ConstantDerivative PlanePosition="Max" Value="10.0"/> -->\n\t\t\t\t\t<!--<Periodic/>-->' \
                      '\n\t\t\t\t</Plane>\n'
            if not flag_2d:
                bc_body += f'\t\t\t\t<Plane Axis="Z">\n\t\t\t\t\t<ConstantValue PlanePosition="Min" Value=' \
                           f'"{item["dirichlet_value"]}"/>\n\t\t\t\t\t<ConstantValue PlanePosition="Max" Value=' \
                           f'"{item["dirichlet_value"]}"/>\n\t\t\t\t\t<!-- Other options are (examples): -->\n\t\t\t\t\t' \
                           f'<!--<ConstantDerivative PlanePosition="Min" Value="10.0"/> -->\n\t\t\t\t\t<!--' \
                           f'<ConstantDerivative PlanePosition="Max" Value="10.0"/> -->\n\t\t\t\t\t<!--<Periodic/>-->' \
                           '\n\t\t\t\t</Plane>\n'
        else:
            bc_body = f'\t\t\t\t<Plane Axis="X">\n\t\t\t\t\t<ConstantDerivative PlanePosition="Min" Value=' \
                      f'"0"/>\n\t\t\t\t\t<ConstantDerivative PlanePosition="Max" Value=' \
                      f'"0"/>\n\t\t\t\t\t<!-- Other options are (examples): -->\n\t\t\t\t\t' \
                      f'<!--<ConstantValue PlanePosition="Min" Value="10.0"/> -->\n\t\t\t\t\t<!--' \
                      f'<ConstantValue PlanePosition="Max" Value="10.0"/> -->\n\t\t\t\t\t<!--<Periodic/>-->' \
                      '\t\t\t\t</Plane>\n' \
                      f'\t\t\t\t<Plane Axis="Y">\n\t\t\t\t\t<ConstantDerivative PlanePosition="Min" Value=' \
                      f'"0"/>\n\t\t\t\t\t<ConstantDerivative PlanePosition="Max" Value=' \
                      f'"0"/>\n\t\t\t\t\t<!-- Other options are (examples): -->\n\t\t\t\t\t' \
                      f'<!--<ConstantValue PlanePosition="Min" Value="10.0"/> -->\n\t\t\t\t\t<!--' \
                      f'<ConstantValue PlanePosition="Max" Value="10.0"/> -->\n\t\t\t\t\t<!--<Periodic/>-->' \
                      '\t\t\t\t</Plane>\n'
            if not flag_2d:
                bc_body += f'\t\t\t\t<Plane Axis="Z">\n\t\t\t\t\t<ConstantDerivative PlanePosition="Min" Value=' \
                           f'"0"/>\n\t\t\t\t\t<ConstantDerivative PlanePosition="Max" Value=' \
                           f'"0"/>\n\t\t\t\t\t<!-- Other options are (examples): -->\n\t\t\t\t\t' \
                           f'<!--<ConstantDerivative PlanePosition="Min" Value="10.0"/> -->\n\t\t\t\t\t<!--' \
                           f'<ConstantDerivative PlanePosition="Max" Value="10.0"/> -->\n\t\t\t\t\t<!--<Periodic/>-->' \
                           '\t\t\t\t</Plane>\n'
        close_bc = "</BoundaryConditions>\n"
        close_field = "</DiffusionField>\n"

        full_field_def = df_str + conc_units + og_D + D_str + og_g + g_str + init_cond_warn + init_cond + het_warning + cells_str + \
                         close_diff_data + bc_head + bc_body + close_bc + close_field
        full_str += full_field_def
    full_str += "</Steppable>\n"
    return full_str


def make_cell_loop(cell_type):
    return f"for cell in self.cell_list_by_type(self.{cell_type.upper()}):\n"


def make_cell_dict(cell_types, secretion_dict):
    for ctype in cell_types:
        loop_start = make_cell_loop(ctype)
        if ctype in secretion_dict.keys():
            type_sec = secretion_dict[ctype]
        else:
            type_sec = None


def convert_secretion_rate(rate, unit, time_conv, pctimeunit, time_convs=_time_convs):
    secretion_comment = ''
    if pctimeunit in unit:  # if it's the same as the "main" time unit
        mcs_rate = rate / time_conv
        return mcs_rate, secretion_comment
    else:
        tu = unit.split('/')[-1]
        if tu not in time_convs.keys():
            message = f"WARNING: Secretion 1/(rate unit) = {tu} not found in {time_convs.keys()}.\nAutomatic conversion" \
                      f" of " \
                      f"this rate is disabled."
            secretion_comment += "\n#" + message.replace("\n", "\n#")
            warnings.warn(message)
            mcs_rate = rate
            return mcs_rate, secretion_comment
        else:
            message = f"WARNING: Secretion 1/(rate unit) = {tu} is not the main PhysiCell time unit {pctimeunit}." \
                      f"\nTherefore, the automatic conversion may be incorrect."
            secretion_comment += "\n#" + message.replace("\n", "\n#")
            warnings.warn(message)
            rate_minutes = rate / time_convs[tu]
            rate_pctime = rate_minutes / time_conv[pctimeunit]
            mcs_rate = rate_pctime / time_conv
            return mcs_rate, secretion_comment

def convert_uptake_rate(rate, unit, time_conv, pctimeunit, time_convs=_time_convs):
    uptake_comment = ''
    if pctimeunit in unit:
        mcs_rate = rate / time_conv
        return mcs_rate, uptake_comment
    else:
        tu = unit.split('/')[-1]
        if tu not in time_convs.keys():
            message = f"WARNING: Uptake 1/(rate unit) = {tu} not found in {time_convs.keys()}.\nAutomatic " \
                      f"conversion of this rate is disabled."
            warnings.warn(message)
            uptake_comment += "#" + message.replace("\n", "\n#")
            mcs_rate = rate
            return mcs_rate, uptake_comment
        else:
            message = f"WARNING: Uptake 1/(rate unit) = {tu} is not the main PhysiCell time unit {pctimeunit}." \
                      f"\nTherefore, the automatic conversion may be incorrect"
            warnings.warn(message)
            uptake_comment += "#" + message.replace("\n", "\n#")

            rate_minutes = rate / time_convs[tu]
            rate_pctime = rate_minutes / time_conv[pctimeunit]
            mcs_rate = rate_pctime / time_conv
            return mcs_rate, uptake_comment


def convert_secretion_data(sec_dict, time_conv, pctimeunit):

    if not sec_dict:
        return {}

    new_sec_dict = sec_dict

    secretion_comment = '#WARNING: PhysiCell has a concept of "target secretion" that CompuCell3D does not. \n#The ' \
                        'translating program attempts to implement it, but it may not be a 1 to 1 conversion.'
    uptake_comment = '#WARNING: To avoid negative concentrations, in CompuCell3D uptake is "bounded." \n# If the amount' \
                     ' that would be uptaken is larger than the value at that pixel,\n# the uptake will be a set ratio ' \
                     'of the amount available.\n# The conversion program uses 1 as the ratio,\n# you may want to ' \
                     'revisit this.'
    for ctype in sec_dict.keys():
        type_sec = sec_dict[ctype]
        new_type_sec = type_sec
        for field, data in type_sec.items():
            unit = data['secretion_unit']

            mcs_secretion_rate, extra_sec_comment = convert_secretion_rate(data['secretion_rate'], unit, time_conv,
                                                                           pctimeunit)
            data['secretion_rate_MCS'] = mcs_secretion_rate
            data['secretion_comment'] = secretion_comment+extra_sec_comment

            mcs_uptake_rate, extra_up_comment = convert_uptake_rate(data['uptake_rate'], data['uptake_unit'], time_conv,
                                                                    pctimeunit)

            data['uptake_rate_MCS'] = mcs_uptake_rate
            data['uptake_comment'] = uptake_comment + extra_up_comment

            new_type_sec[field] = data
        new_sec_dict[ctype] = new_type_sec
    return new_sec_dict



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

    example_path = r"./example_pcxml/" + "annotated_cancer_immune3D_flat.xml"

    print(f"Loading {example_path}")
    with open(example_path, 'r') as f:
        xml_raw = f.read()
    pcdict = x2d.parse(xml_raw)['PhysiCell_settings']

    print("Generating <Metadata/>")
    metadata_str, n_threads = make_metadata(pcdict)

    print("Generating <Potts/>")
    potts_str, pcdims, ccdims, pctime, cctime = make_potts(pcdict)

    print("Generating <Plugin CellType/>")
    ct_str, wall, cell_types, = make_cell_type_plugin(pcdict)

    constraints = get_cell_constraints(pcdict, ccdims[4], cctime[2])

    with open(os.path.join(out_sim_f, "extra_definitions.py"), 'w+') as f:
        f.write("cell_constraints=" + str(constraints) + "\n")

    print("Generating <Plugin Contact/>")
    contact_plug = make_contact_plugin(cell_types)

    extra = extra_for_testing(cell_types, ccdims[0], ccdims[1], ccdims[2])

    print("parsing micro environment")
    d_elements = get_microenvironment(pcdict, ccdims[4], pcdims[3], cctime[2], pctime[1])

    print("Generating diffusion plugin")
    diffusion_string = make_diffusion_plug(d_elements, cell_types, False)

    secretion_dict = get_secretion(pcdict)

    conv_sec = convert_secretion_data(secretion_dict, cctime[2], pctime[1]) 

    constraint_step = steppable_gen.generate_constraint_steppable(cell_types,
                                                                  [constraints,
                                                                   conv_sec])

    print("Merging")
    cc3dml = "<CompuCell3D>\n"
    cc3dml += metadata_str + potts_str + ct_str + contact_plug + diffusion_string + '\n' + extra + \
              "\n</CompuCell3D>\n"

    print(f"Creating {out_sim_f}/test.xml")
    with open(os.path.join(out_sim_f, "test.xml"), "w+") as f:
        f.write(cc3dml)

    print("Copying python files")

    sh.copy(r'./base_cc3d_python_scripts/test.py', out_sim_f)
    sh.copy(r'./base_cc3d_python_scripts/testSteppables.py', out_sim_f)

    print("______________\nDONE!!")
    # print(cc3dml)
