# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 09:45:52 2022

@author: Juliano Ferrari Gianlupi
"""


# import sys
# import string
# import copy
import os
import warnings
import shutil as sh
import xmltodict as x2d

from itertools import combinations

time_convs = {}


def get_dims(pcdict):
    xmin = float(pcdict['domain']['x_min']) if "x_min" in pcdict['domain'].keys() else None
    xmax = float(pcdict['domain']['x_max']) if "x_max" in pcdict['domain'].keys() else None

    ymin = float(pcdict['domain']['y_min']) if "y_min" in pcdict['domain'].keys() else None
    ymax = float(pcdict['domain']['y_max']) if "y_max" in pcdict['domain'].keys() else None

    zmin = float(pcdict['domain']['z_min']) if "z_min" in pcdict['domain'].keys() else None
    zmax = float(pcdict['domain']['z_max']) if "z_max" in pcdict['domain'].keys() else None

    units = pcdict['overall']['space_units'] if 'overall' in pcdict.keys() and \
                                                'space_units' in pcdict['overall'].keys() else 'micron'

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

    cc3dspaceunitstr = f"1 pixel = {cc3dds} {units}"

    return ((xmin, xmax), (ymin, ymax), (zmin, zmax), units), \
           (cc3dx, cc3dy, cc3dz, cc3dspaceunitstr, cc3dds)


def get_time(pcdict):
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

    mechdt = float(pcdict['overall']['dt_mechanics']['#text']) if "dt_mechanics" in pcdict['overall'].keys() else None

    steps = round(mt / mechdt)

    cc3ddt = 1 / (mt / steps)  # MCS/unit

    cc3dtimeunitstr = f"1 MCS = {cc3ddt} {time_unit}"

    # timeconvfact = 1/cc3ddt

    return (mt, time_unit, mechdt), (steps, cc3dtimeunitstr, cc3ddt)


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

    # sys.exit()

    with open(os.path.join(out_sim_f, "extra_definitions.py"), 'w+') as f:
        f.write("cell_constraints=" + str(constraints) + "\n")

    print("Generating <Plugin Contact/>")
    contact_plug = make_contact_plugin(cell_types)

    extra = extra_for_testing(cell_types, ccdims[0], ccdims[1], ccdims[2])

    print("Merging")
    cc3dml = "<CompuCell3D>\n"
    cc3dml += metadata_str + potts_str + ct_str + contact_plug + '\n' + extra + \
              "\n</CompuCell3D>\n"

    print(f"Creating {out_sim_f}/test.xml")
    with open(os.path.join(out_sim_f, "test.xml"), "w+") as f:
        f.write(cc3dml)

    print("Copying python files")

    sh.copy(r'./base_cc3d_python_scripts/test.py', out_sim_f)
    sh.copy(r'./base_cc3d_python_scripts/testSteppables.py', out_sim_f)

    print("______________\nDONE!!")
    # print(cc3dml)

