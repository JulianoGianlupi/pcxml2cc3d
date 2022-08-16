# -*- coding: utf-8 -*-
"""
Created on Mon Aug 15 09:45:52 2022

@author: Juliano Ferrari Gianlupi
"""


import xml.etree.ElementTree as ET
import sys
import string
import copy
import os
import warnings
import shutil as sh

from itertools import combinations

time_convs = {}

def get_tags(root):
    return [c.tag for c in root.iter()]

def get_dims(tags, root):
    xml_root = root #todo: clean
    
    # 
    
    xmin = float(next(xml_root.iter("x_min")).text) if "x_min" in tags \
        else None
    xmax = float(next(xml_root.iter("x_max")).text) if "x_max" in tags \
        else None
    
    ymin = float(next(xml_root.iter("y_min")).text) if "y_min" in tags \
        else None
    ymax = float(next(xml_root.iter("y_max")).text) if "y_max" in tags \
        else None
    
    zmin = float(next(xml_root.iter("z_min")).text) if "x_min" in tags \
        else None
    zmax = float(next(xml_root.iter("z_max")).text) if "z_max" in tags \
        else None
    
    units = "micron" if "space_units" not in tags else \
        next(xml_root.iter("space_units")).text
        
    # the dx/dy/dz tags mean that for every voxel there are dx space-units.
    # therefore [dx] = [space-unit/voxel]. Source: John Metzcar
    
    dx = float(next(xml_root.iter("dx")).text) if "dx" in tags else 1
    dy = float(next(xml_root.iter("dx")).text) if "dy" in tags else 1
    dz = float(next(xml_root.iter("dx")).text) if "dz" in tags else 1
    print(dx,dy,dz, type(dx),type(dy),type(dz),)
    if not dx == dy == dz:
        message="WARNING! Physicell's dx/dy/dz are not all the same: "\
            f"dx={dx}, dy={dy}, dz={dz}\n"\
                f"Using {max(min(dx,dy,dz), 1)}"
        warnings.warn(message)
    
    ds = max(min(dx,dy,dz), 1)
    
    diffx = 1 if xmin is None or xmax is None else round(xmax-xmin)
    diffy = 1 if ymin is None or ymax is None else round(ymax-ymin)
    diffz = 1 if zmin is None or zmax is None else round(ymax-zmin)
    
    cc3dx = round(diffx/ds) 
    cc3dy = round(diffy/ds)
    cc3dz = round(diffz/ds)
    
    cc3dds = cc3dx/diffx # pixel/unit
    
    cc3dspaceunitstr = f"1 pixel = {cc3dds} {units}"
    
    return ((xmin, xmax), (ymin,ymax), (zmin,zmax), units),\
            (cc3dx, cc3dy, cc3dz, cc3dspaceunitstr, cc3dds)

def get_time(tags, root):
    xml_root = root # todo: clean
    mt = float(next(xml_root.iter("max_time")).text) if "max_time" in tags \
        else 100000
    mtunit = next(xml_root.iter("max_time")).attrib["units"] if "max_time" in tags \
        else None
        
    time_unit =  next(xml_root.iter("time_units")).text if "time_units" in tags \
        else None
        
    if mtunit != time_unit:
        message = f"Warning: Psysicell time units in "\
                "\n`<overall>\n\t<max_time units=...`\ndiffers from\n"\
                    f"`<time_units>unit</time_units>`.\nUsing: {time_unit}"
        warnings.warn(message)
    mechdt = float(next(xml_root.iter("dt_mechanics")).text) if "dt_mechanics"\
        in tags else 1
    
    steps = round(mt/mechdt)
    
    cc3ddt = 1/(mt/steps) # MCS/unit
    
    cc3dtimeunitstr = f"1 MCS = {cc3ddt} {time_unit}"
    
    # timeconvfact = 1/cc3ddt
    
    return (mt, time_unit, mechdt), (steps, cc3dtimeunitstr, cc3ddt)

def get_parallel(tags, root):
    return int(next(root.iter("omp_num_threads")).text) if "dt_mechanics"\
        in tags else 1

def make_potts(tags, root):
    '''
    

    Parameters
    ----------
    tags : List
        A list of all xml tags in the PhysiCell xml file.
    root : xml.etree.ElementTree.Element
        the xml that has g.

    Returns
    -------
    potts_str: string
        the generated CC3D XML Potts block.
    pcdim: tupple
        the extracted spatial dimensions from PhysiCell XML
    
    ccdims: tupple
    
    pctime: tupple
    
    cctime: tupple

    '''
    
    
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
    pcdims, ccdims = get_dims(tags, root) # ((xmin, xmax), (ymin,ymax), (zmin,zmax), units), (cc3dx, cc3dy, cc3dz, cc3dspaceunitstr, cc3dds)
    
    
    
    # space_units_str = f'"1 pixel = 1 {space_units}"'
    
    pctime, cctime = get_time(tags, root)
    
    
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

def make_metadata(tags, root, out=100):
    
    threads = get_parallel(tags, root)
    
    metadata = f'''
<Metadata>
  <!-- Basic properties simulation -->
  <NumberOfProcessors>{threads}</NumberOfProcessors>
  <DebugOutputFrequency>{out}</DebugOutputFrequency>
  <!-- <NonParallelModule Name="Potts"/> -->
</Metadata>\n'''
    
    return metadata, threads

def get_boundary_wall(tags, root):
    
    if "virtual_wall_at_domain_edge" in tags:
        wall = next(root.iter("virtual_wall_at_domain_edge")).text.upper()
        if wall == "TRUE":
            return True 
        else:
            return False
    else:
        return False
    return False

def make_cell_type_tags(tags,root):
    
    s = ''
    cell_types = []
    idx = 1
    for child in root.iter("cell_definition"):
        # print(child.tag, child.attrib, child.text)
        name = child.attrib['name'].replace(" ", "_")
        cell_types.append(name)
        ctt = f'\t<CellType TypeId="{idx}" TypeName="{name}"/>\n'
        s += ctt
        idx+=1
        
    create_wall = get_boundary_wall(tags, root)
    
    if create_wall:
        s += f'\t<CellType Freeze="" TypeId="{idx}" TypeName="WALL"/>\n'
        cell_types.append("WALL")
    
    return s, create_wall, cell_types

def make_cell_type_plugin(tags,root):
    
    ct_str = '\n<Plugin Name="CellType">\n\t'\
        '<CellType TypeId="0" TypeName="Medium"/>\n'
    typesstr, wall, cell_types = make_cell_type_tags(tags, root)
    
    ct_str += typesstr
    ct_str += '</Plugin>'
    
    return ct_str, wall, cell_types

def make_cc3d_file(name=None):
    
    if name is None:
        cc3d = '''
<Simulation version="4.3.0">
   <XMLScript Type="XMLScript">Simulation/test.xml</XMLScript>
    <PythonScript Type="PythonScript">Simulation/test.py</PythonScript> 
    <Resource Type="Python">Simulation/testSteppables.py</Resource> 
</Simulation>\n'''
        return cc3d
    else:
        cc3d = '''
<Simulation version="4.3.0">
   <XMLScript Type="XMLScript">Simulation/{name}.xml</XMLScript>
   <PythonScript Type="PythonScript">Simulation/{name}.py</PythonScript>
   <Resource Type="Python">Simulation/{name}Steppables.py</Resource>
</Simulation>\n'''
        return cc3d

def make_contact_plugin(celltypes):
    # replace ASAP with contact flex
    
    
    
    combs = list(combinations(celltypes, 2))
    
    for t in celltypes:
        combs.append((t,t))
    
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
    
    beg='''<Steppable Type="UniformInitializer">
\t<!-- Initial layout of cells in the form of rectangular slab -->
\t<Region>\n'''
    
    box_min = f'\t\t<BoxMin x="{max(1, xmax-50)}" y="{max(1, ymax-50)}" z="{max(1, xmax-50)}"/>\n'
    box_max = f'\t\t<BoxMax x="{xmax-10}" y="{ymax-10}" z="{zmax-10}"/>\n'
    
    gap = "\t\t<Gap>0</Gap>\n\t\t<Width>7</Width>\n"

    types = ''
    for t in celltypes: 
        if t.upper()=="WALL":
            continue
        types += f"{t},"
    
    types = types[:-1]
    
    types = "\t\t<Types>"+types+"</Types>\n"
    
    end='''\t</Region>
</Steppable>'''
    
    return beg+box_min+box_max+gap+types+end
    

if __name__=="__main__":
    
    
    
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
    
    example_path = r"./example_pcxml/"+"annotated_cancer_immune3D_flat.xml"
    
    print(f"Loading {example_path}")
    tree = ET.parse(example_path)
    xml_root = tree.getroot()
    
    
    
    print("Getting PhysiCell XML tags")
    tags = get_tags(xml_root)
    
    print("Generating <Metadata/>")
    metadata_str, n_threads = make_metadata(tags, xml_root)
    
    print("Generating <Potts/>")
    potts_str, pcdims, ccdims, pctime, cctime = make_potts(tags, xml_root)
    
    
    print("Generating <Plugin CellType/>")
    ct_str, wall, cell_types = make_cell_type_plugin(tags, xml_root)
    
    print("Generating <Plugin Contact/>")
    contact_plug = make_contact_plugin(cell_types)
    
    extra = extra_for_testing(cell_types, ccdims[0], ccdims[1], ccdims[2])
    
    print("Merging")
    cc3dml = "<CompuCell3D>\n"
    cc3dml += metadata_str+potts_str+ct_str+contact_plug +'\n'+extra+\
        "\n</CompuCell3D>\n"
    
    print(f"Creating {out_sim_f}/test.xml")
    with open(os.path.join(out_sim_f, "test.xml"), "w+") as f:
        f.write(cc3dml)
        
    print("Copying python files")
    
    sh.copy(r'./base_cc3d_python_scripts/test.py', out_sim_f)
    sh.copy(r'./base_cc3d_python_scripts/testSteppables.py', out_sim_f)
    
    print("______________\nDONE!!")
    # print(cc3dml)
    
    for child in xml_root.iter():
        break
        print(child.tag, child.attrib, child.text)
        if child.tag == "max_time":
            
            break
    