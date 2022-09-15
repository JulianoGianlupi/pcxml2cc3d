from itertools import combinations

from cc3d_xml_gen.get_physicell_data import get_dims, get_time, get_parallel, get_boundary_wall

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

def make_potts(pcdict):
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


def make_cell_type_plugin(pcdict):
    ct_str = '\n<Plugin Name="CellType">\n\t' \
             '<CellType TypeId="0" TypeName="Medium"/>\n'
    typesstr, wall, cell_types = make_cell_type_tags(pcdict)

    ct_str += typesstr
    ct_str += '</Plugin>'

    return ct_str, wall, cell_types


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