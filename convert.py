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

time_convs = {}

def get_tags(root):
    return [c.tag for c in root.iter()]

def get_dims(tags, root):
    xml_root = root #todo: clean
    
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
    
    cc3dx = 1 if xmin is None or xmax is None else str(round(xmax-xmin))
    cc3dy = 1 if ymin is None or ymax is None else str(round(ymax-ymin))
    cc3dz = 1 if zmin is None or zmax is None else str(round(ymax-zmin))
    
    return ((xmin, xmax), (ymin,ymax), (zmin,zmax)),\
            (cc3dx, cc3dy, cc3dz)

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




if __name__=="__main__":
    print("Running test")
    example_path = r"./example_pcxml/"+"cancer_immune3D_flat.xml"
    
    tree = ET.parse(example_path)
    xml_root = tree.getroot()
    
    tags = get_tags(xml_root)
    
    
    pcdims, ccdims = get_dims(tags, xml_root)
    
    pctime, cctime = get_time(tags, xml_root)
    
    
    potts_str = f""" # still need to implement space units
    <Potts>
      <!-- Basic properties of CPM (GGH) algorithm -->
      <Dimensions x="{ccdims[0]}" y="{ccdims[1]}" z="{ccdims[2]}"/>
      <Time_Units>"{cctime[1]}"</Time_Units>
      <MCS_to_time units="{pctime[1]}/MCS">{cctime[2]}</MCS_to_time>
      <Steps>{cctime[0]}</Steps>
      <!-- As the frameworks of CC3D and PhysiCell are very different -->
      <!-- PC doesn't have some concepts that CC3D does. Temperature is one of -->
      <!-- them, so the translation script leaves its tunning as an exercise-->
      <!-- for the reader -->
      <Temperature>10.0</Temperature>
      <!-- same deal for neighbor order -->
      <NeighborOrder>1</NeighborOrder>
      <!-- <Boundary_x>Periodic</Boundary_x> -->
      <!-- <Boundary_y>Periodic</Boundary_y> -->
   </Potts>
    """
    print(potts_str)
    
    for child in xml_root.iter():
        break
        print(child.tag, child.attrib, child.text)
        if child.tag == "max_time":
            
            break
    