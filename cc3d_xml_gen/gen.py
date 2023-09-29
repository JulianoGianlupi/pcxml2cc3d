import warnings
from itertools import combinations
from math import ceil

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


def make_potts(pcdims, ccdims, pctime, cctime):
    """
    Generate a Potts CC3D XML string with the given parameters.

    Parameters
    ----------
    pcdims : tuple
        Tuple of integers representing the dimensions of space in PhysiCell
    ccdims : tuple
        Tuple of integers representing the dimensions of space in CC3D
    pctime : tuple
        Tuple with the time unit from PhysiCell
    cctime : tuple
        Tuple with the time parameters for CC3D

    Returns
    -------
    str
        Potts XML string with the given parameters.

    """
    potts_str = f""" 
<Potts>
   <!-- Basic properties of CPM (GGH) algorithm -->
   <Space_Units>{ccdims[3]}</Space_Units>
   <Pixel_to_Space units="pixel/{pcdims[3]}" id = "pixel_to_space">{ccdims[4]}</Pixel_to_Space>
   <Dimensions x="{ccdims[0]}" y="{ccdims[1]}" z="{ccdims[2]}"/>
   <Time_Units>{cctime[1]}</Time_Units>
   <MCS_to_Time units="MCS/{pctime[1]}" id = "mcs_to_time">{cctime[2]}</MCS_to_Time>
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

    return potts_str


def make_metadata(pcdict, out=100):
    """Generates the metadata CC3D XML block"""
    threads = get_parallel(pcdict)

    metadata = f'''
<Metadata>
  <!-- Basic properties simulation -->
  <NumberOfProcessors>{1}</NumberOfProcessors>
  <DebugOutputFrequency>{out}</DebugOutputFrequency>
  <!-- <NonParallelModule Name="Potts"/> -->
</Metadata>\n'''

    return metadata, threads


def make_cell_type_plugin(pcdict):
    """
    Makes the cell type plugin for CC3D

    Passes pcdict to make_cell_type_tags to generate the cell types and returns the results

    :param pcdict: Dictionary created from parsing PhysiCell XML
    :return ct_str, wall, cell_types: string setting the cell type plugin for cc3d's XML, bool for the presence of a
    wall cell type, list of cell types
    """
    ct_str = '\n<Plugin Name="CellType">\n\t' \
             '<CellType TypeId="0" TypeName="Medium"/>\n'
    typesstr, wall, cell_types = make_cell_type_tags(pcdict)

    ct_str += typesstr
    ct_str += '</Plugin>'

    return ct_str, wall, cell_types


def make_cell_type_tags(pcdict):
    """
    Parses the PhysiCell dictionary to fetch the cell type names, generates the internal part of the cell type plugin


    :param pcdict: Dictionary created from parsing PhysiCell XML
    :return s, create_wall, cell_types: string for the cell type plugin, bool for the existance of a wall cell type,
    list of cell type names
    """
    s = ''
    cell_types = []
    idx = 1

    create_wall = get_boundary_wall(pcdict)

    if create_wall:
        s += f'\t<CellType Freeze="" TypeId="{idx}" TypeName="WALL"/>\n'
        cell_types.append("WALL")
        idx += 1

    if 'cell_definitions' not in pcdict.keys():
        name = "CELL"
        cell_types.append(name)
        ctt = f'\t<CellType TypeId="{idx}" TypeName="{name}"/>\n'
        s += ctt
        return s, create_wall, cell_types

    for child in pcdict['cell_definitions']['cell_definition']:
        # print(child.tag, child.attrib, child.text)
        name = child['@name'].replace(" ", "_")
        cell_types.append(name)
        ctt = f'\t<CellType TypeId="{idx}" TypeName="{name}"/>\n'
        s += ctt

        idx += 1

    return s, create_wall, cell_types


def make_cc3d_file(name=None):
    """
    Generates the string for the .cc3d file and returns the names of the simulation files

    Based on the name parameter, this function generates the XML for the .cc3d file and the file
    names for all the simulation files

    :param name: string, simulation name
    :return: string for .cc3d file, name of the xml file, name of the main python file, name of the steppables file
    """
    if name is None:
        xml_name = "test.xml"
        main_py_name = "main_test.py"
        steppables_py_name = "steppable_test.py"
    else:
        xml_name = f"{name}.xml"
        main_py_name = f"{name}.py"
        steppables_py_name = f"{name}Steppables.py"
    cc3d = f'''
<Simulation version="4.3.0">
   <XMLScript Type="XMLScript">Simulation/{xml_name}</XMLScript>
   <PythonScript Type="PythonScript">Simulation/{main_py_name}</PythonScript>
   <Resource Type="Python">Simulation/{steppables_py_name}</Resource>
   <Resource Type="Python">Simulation/extra_definitions.py</Resource> 
</Simulation>\n'''
    return cc3d, xml_name, main_py_name, steppables_py_name


def make_contact_plugin(celltypes):
    """
    Create a contact plugin configuration XML for CompuCell3D simulations.

    The function generates the configuration XML for the contact plugin used in a CompuCell3D simulation. The plugin
    defines the energy associated with different types of cell-cell interactions (defaulted to 10) and with the
    interaction between cells and the surrounding medium (defaulted to 10). The configuration XML is returned as a
    string.

    The function first creates all combinations of the given cell types using the combinations function from the
    itertools module, and then appends the self-interactions (cell-cell interactions of the same type) to the list of
    combinations. The list is then reversed, so that the self-interactions appear first.

    The function then creates the configuration XML as a string, with 10 as the placeholder energy. The NeighborOrder
    parameter is set to 3.

    Parameters:
    ----------
        celltypes : list
            A list of strings representing the cell types in the simulation.

    Returns:
    -------
        contact_plug : str
            A string containing the configuration for the contact plugin.
    """
    # todo: replace ASAP with contact flex

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
        me += f'\t<Energy Type1="Medium" Type2="{t}">5.0</Energy>\n'

    # 2 make the combination energies

    ce = ""
    for t1, t2 in combs:
        ce += f'\t<Energy Type1="{t1}" Type2="{t2}">5.0</Energy>\n'

    contact_plug += me + ce + "\t<NeighborOrder>3</NeighborOrder>\n</Plugin>"
    return contact_plug


def make_diffusion_FE(diffusing_elements, celltypes, flag_2d):
    """
    Converts a dictionary of diffusion properties into a CC3D DiffusionSolverFE XML configuration string. T

    This function generates an XML string that can be used to configure CC3D's DiffusionSolverFE. It takes three
    arguments: diffusing_elements, celltypes, and flag_2d. diffusing_elements is a dictionary of the diffusing
    elements, where each key is the name of the diffusing element and the corresponding value is another dictionary
    containing the properties of that element, such as the diffusion constant and initial concentration. celltypes is a
    list of the cell types, and flag_2d is a boolean indicating whether the simulation is in two dimensions or not.

    The function loops through each diffusing element in the dictionary and generates a string with information about
    the diffusion field, including its name, diffusion data (such as diffusion and decay constants), initial
    concentration, and boundary conditions. It then concatenates these strings together to create the full XML string.

    Parameters
    ----------
    diffusing_elements : dict
        A dictionary of diffusion properties. Each key represents a diffusing element and contains a nested dictionary
        with keys "D", "gamma", "concentration_units", "D_w_units", "D_og_unit", "gamma_w_units", "gamma_og_unit",
        "use_steady_state", "initial_condition", "dirichlet", "dirichlet_value"
    celltypes : list
        A list of cell types.
    flag_2d : bool
        A boolean value indicating whether the simulation is in 2D or 3D.

    Returns
    -------
    str
        A string representation of the DiffusionSolverFE XML configuration file.

    """
    header = f'\n\n\t<Steppable Type="DiffusionSolverFE">\n\t\t<!-- The conversion uses DiffusionSolverFE and' \
             f' SteadyStateDiffusionSolver ' \
             f'by default. You may ' \
             f'wish to use another diffusion solver-->\n'

    full_str = header

    for key, item in diffusing_elements.items():

        if item["use_steady_state"]:
            continue

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
        if item['dirichlet'].upper() != "FALSE":
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
                           f'<ConstantDerivative PlanePosition="Max" Value="0.0"/> -->\n\t\t\t\t\t<!--<Periodic/>-->' \
                           '\t\t\t\t</Plane>\n'
        close_bc = "</BoundaryConditions>\n"
        close_field = "</DiffusionField>\n"

        full_field_def = df_str + conc_units + og_D + D_str + og_g + g_str + init_cond_warn + init_cond + het_warning + cells_str + \
                         close_diff_data + bc_head + bc_body + close_bc + close_field
        full_str += full_field_def
    full_str += "</Steppable>\n"
    return full_str


def make_diffusion_steady(diffusing_elements, flag_2d):
    """
    Creates a steady-state diffusion solver configuration for CC3D simulations based on the given diffusing elements.

    This function generates an XML string that can be used to configure CC3D's steady state diffusion solver. It takes
    two arguments: diffusing_elements and flag_2d. diffusing_elements is a dictionary of the diffusing
    elements, where each key is the name of the diffusing element and the corresponding value is another dictionary
    containing the properties of that element, such as the diffusion constant and initial concentration. flag_2d is a
    boolean indicating whether the simulation is in two dimensions or not.

    The function loops through each diffusing element in the dictionary and generates a string with information about
    the diffusion field, including its name, diffusion data (such as diffusion and decay constants), initial
    concentration, and boundary conditions. It then concatenates these strings together to create the full XML string.

    Parameters
    ----------
    diffusing_elements : dict
        diffusing_elements : dict
        A dictionary of diffusion properties. Each key represents a diffusing element and contains a nested dictionary
        with keys "D", "gamma", "concentration_units", "D_w_units", "D_og_unit", "gamma_w_units", "gamma_og_unit",
        "use_steady_state", "initial_condition", "dirichlet", "dirichlet_value"

    flag_2d : bool
        A boolean indicating whether to use 2D or 3D solver.

    Returns
    -------
    str
        A string containing the configuration for the steady-state diffusion solver in CC3D simulations.
    """
    if flag_2d:
        header = f'\n\n\t<Steppable Type="SteadyStateDiffusionSolver2D">\n\t\t<!-- The conversion uses ' \
                 f'DiffusionSolverFE and' \
                 f' SteadyStateDiffusionSolver ' \
                 f'by default. You may ' \
                 f'wish to use another diffusion solver-->\n'
    else:
        header = f'\n\n\t<Steppable Type="SteadyStateDiffusionSolver">\n\t\t<!-- The conversion uses ' \
                 f'DiffusionSolverFE and' \
                 f' SteadyStateDiffusionSolver ' \
                 f'by default. You may ' \
                 f'wish to use another diffusion solver-->\n'

    full_str = header

    for key, item in diffusing_elements.items():
        if not item["use_steady_state"]:
            continue

        name = key.replace(" ", "_")

        df_str = f'\t\t<DiffusionField Name="{name}">\n\t\t\t<DiffusionData>\n\t\t\t\t<FieldName>{name}</FieldName>\n'
        conc_units = f'\t\t\t\t<Concentration_units>{item["concentration_units"]}</Concentration_units>\n'
        og_D = f'\t\t\t\t<Original_diffusion_constant D="{item["D_w_units"]}" units= "{item["D_og_unit"]}"/>\n'

        D_str = f'\t\t\t\t<DiffusionConstant>{item["D"]}</DiffusionConstant>\n'
        og_g = f'\t\t\t\t<Original_decay_constant gamma="{item["gamma_w_units"]}" units= "{item["gamma_og_unit"]}"/>\n'
        g_str = f'\t\t\t\t<DecayConstant>{item["gamma"]}</DecayConstant>\n'

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

        close_diff_data = "\t\t\t</DiffusionData>\n"

        # boundary conditions

        bc_head = '\t\t\t<BoundaryConditions>\n\t\t\t\t<!-- PhysiCell has either Dirichlet boundary conditions (i.e. ' \
                  'constant ' \
                  'value) -->\n\t\t\t\t<!-- or "free floating" boundary conditions (i.e., constant flux = 0). -->' \
                  '\n\t\t\t\t<!-- CC3D ' \
                  'allows ' \
                  'for more control of boundary conditions, you may want to revisit the issue. -->\n'
        if item['dirichlet'].upper() != "FALSE":
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

        full_field_def = df_str + conc_units + og_D + D_str + og_g + g_str + init_cond_warn + init_cond + \
                         close_diff_data + bc_head + bc_body + close_bc + close_field
        full_str += full_field_def

    full_str += "</Steppable>\n"
    return full_str


def determine_diffusion_existence(diffusing_elements):
    steadys = []
    for _, item in diffusing_elements.items():
        steadys.append(item["use_steady_state"])
    if not bool(len(steadys)):
        return False, False
    steady = any(steadys)
    regular = any([not el for el in steadys])
    return regular, steady


def make_diffusion_plug(diffusing_elements, celltypes, flag_2d):
    """
    Generates the XML for the diffusion steppables in CC3D a diffusion plug, combining the finite element (FE) solver
    and the steady-state solver.

    Parameters
    ----------
    diffusing_elements : dict
        Dictionary of diffusing elements and their parameters
        Each key in the dictionary represents a diffusing element and its value is a dictionary with the following keys:
        - use_steady_state : bool
            Whether or not to use steady-state diffusion solver for this element.
        - concentration_units : str
            The concentration units of this element.
        - D_w_units : float
            The diffusion constant with units of concentration^2/time.
        - D_og_unit : str
            The original unit of the diffusion constant.
        - D : float
            The diffusion constant without units.
        - gamma_w_units : float
            The decay constant with units of 1/time.
        - gamma_og_unit : str
            The original unit of the decay constant.
        - gamma : float
            The decay constant without units.
        - initial_condition : str
            The initial concentration expression for this element.
        - dirichlet : bool
            Whether or not to use Dirichlet boundary conditions for this element.
        - dirichlet_value : float
            The value of Dirichlet boundary condition for this element.
    celltypes : list
        List of cell types
    flag_2d : bool
        Whether the simulation is in 2D

    Returns
    -------
    str
        The combined string of the FE and steady-state solvers
    """

    use_regular, use_steady = determine_diffusion_existence(diffusing_elements)

    if use_regular:
        FE_solver = make_diffusion_FE(diffusing_elements, celltypes, flag_2d)
    else:
        FE_solver = ""
    if use_steady:
        steady_state_solver = make_diffusion_steady(diffusing_elements, flag_2d)
    else:
        steady_state_solver = ""

    return FE_solver + steady_state_solver


def get_chemotatic_fields(taxis_dict):
    taxis_fields = []
    for ctype, taxis_data in taxis_dict.items():
        for key in taxis_data.keys():
            if key != 'chemotaxis_dict':
                taxis_fields.append(key)
    return list(set(taxis_fields))


def make_chemotaxis(taxis_dict):
    if not bool(taxis_dict):
        return ""
    else:
        taxis_fields = get_chemotatic_fields(taxis_dict)
        plug = '<Plugin Name="Chemotaxis">\n'
        for name in taxis_fields:
            plug += f'<ChemicalField Name="{name}"/>\n'
        plug += '</Plugin>'
        return plug


def make_secretion(secretion_dict):
    if bool(secretion_dict):
        return '\n<Plugin Name="Secretion"/>\n'
    else:
        return ""


def make_volume():
    return '\n<Plugin Name="Volume"/>\n'


def make_cell_loop(cell_type):
    return f"for cell in self.cell_list_by_type(self.{cell_type.upper()}):\n"


def make_cell_dict(cell_types, secretion_dict):
    for ctype in cell_types:
        loop_start = make_cell_loop(ctype)
        if ctype in secretion_dict.keys():
            type_sec = secretion_dict[ctype]
        else:
            type_sec = None


def reconvert_cc3d_dims(ccdims, ratio, is_2D):
    """
    Convert the number of pixels in each dimension of a 3D domain based on a given ratio.

    Parameters:
    -----------
        ccdims : tuple
            A tuple containing the number of pixels in each dimension of a 3D domain, the pixel --
            real unit relationship,
            and the pixel -- real unit ratio.

        ratio : int
            The ratio to multiply the number of pixels in each dimension by.

    Returns:
    --------
        tuple: A new tuple containing the updated number of pixels in each dimension,
              the updated pixel -- real unit relationship, and the updated pixel -- real unit ratio.
              :param is_2D:
              :type is_2D:
    """
    # ccdims is a tupple of: int x pixels, int y pixels, int z pixels,
    # string showing the pixel -- real unit relationship, the pixel -- real unit ratio

    ccdims = list(ccdims)

    if not is_2D:
        number_pixels = [ratio * ccdims[0], ratio * ccdims[1], ratio * ccdims[2]]
    else:
        number_pixels = [ratio * ccdims[0], ratio * ccdims[1], 1]

    old_pixel_unit_ratio = ccdims[4]

    # pixel` = ratio*pixel = ratio * conv * unit
    new_pixel_unit_ratio = ratio * old_pixel_unit_ratio

    new_string = ccdims[3].replace(str(old_pixel_unit_ratio), str(new_pixel_unit_ratio))
    new_ccdims = (number_pixels[0], number_pixels[1], number_pixels[2], new_string, new_pixel_unit_ratio, ccdims[-1],
                  is_2D)
    return new_ccdims


def reconvert_cell_volume_constraints(con_dict, ratio, minimum_volume):
    """
    Parameters
    ----------
    con_dict : dict
       Dictionary of constraints to be converted.
       Each key in the dictionary corresponds to a type of constraint,
       and its value is a dictionary with keys like 'volume', 'surface', etc.
    ratio : int
       The ratio between the old voxel size and the new voxel size.
       All volume constraints are multiplied by this factor.
    minimum_volume : int
       The minimum value to replace any None values in the 'volume' key.

    Returns
    -------
    dict
       A new dictionary of constraints, where the 'volume' key values
       have been converted to the new voxel size, and any None values
       have been replaced with the `minimum_volume`.
    """
    new_con = {}
    for ctype, const in con_dict.items():
        new_con[ctype] = const
        if const['volume']["volume (pixels)"] is None:
            new_con[ctype]['volume']["volume (pixels)"] = minimum_volume
        else:
            new_con[ctype]['volume']["volume (pixels)"] = ratio * const['volume']["volume (pixels)"]
    return new_con


def reconvert_spatial_parameters_with_minimum_cell_volume(constraints, ccdims, pixel_volumes, minimum_volume):
    """
    Convert spatial parameters `ccdims` and `constraints` based on a minimum cell volume.

    The reconvert_spatial_parameters_with_minimum_cell_volume function takes in four arguments: constraints, ccdims,
    pixel_volumes, and minimum_volume. It returns a tuple containing the converted ccdims and constraints.

    The function first filters out any None values from the pixel_volumes list and determines the minimum volume from
    the remaining values. It then calculates a reconvert_ratio by dividing the minimum_volume by the
    minimum_converted_volume and taking the ceiling of the result.

    The ccdims argument is then converted using the reconvert_cc3d_dims function with the reconvert_ratio. The
    constraints argument is also converted using the reconvert_cell_volume_constraints function with the
    reconvert_ratio and minimum_volume as arguments.

    The converted ccdims and constraints are then returned as a tuple.

    Parameters:
    -----------
        constraints : dict
            The previously converted dictionary of cell constraints
        ccdims : tupple
            A tuple of the previously converted cc3d space parameters
        pixel_volumes : list
            A list of volumes of cells in pixels
        minimum_volume : int
            The minimum volume required for the cell in pixels

    Returns:
        - Tuple[Tuple, Dict]: A tuple containing the converted `ccdims` and `constraints`.
    """
    is_2D = ccdims[6]
    px_vols = [px for px in pixel_volumes if px is not None]
    minimum_converted_volume = min(px_vols)

    reconvert_ratio = ceil(minimum_volume / minimum_converted_volume)

    ccdims = reconvert_cc3d_dims(ccdims, reconvert_ratio, is_2D)

    constraints = reconvert_cell_volume_constraints(constraints, reconvert_ratio, minimum_volume)

    return ccdims, constraints


def decrease_domain(ccdims, max_volume=150 ** 3):
    """
    Decrease the size of a 3D domain if its volume exceeds a maximum value.

    This function calculates the volume of the domain based on its current dimensions,
    and compares it to a maximum volume value. If the volume exceeds the maximum, the
    function decreases the size of the domain by scaling the dimensions proportionally
    to the original dimensions. If the original dimensions are not uniform, the scaling
    is applied to the default side length (the cubic root of the maximum volume) times
    the proportions of the original dimensions. If the original dimensions are uniform,
    the scaling is applied to all dimensions uniformly. The new dimensions are returned
    along with a boolean indicating whether the domain was truncated. If the domain was
    truncated, a warning message is issued, indicating that the initial conditions defined
    in PhysiCell may be affected.

    Parameters:
    -----------
        ccdims : tuple
            tuple of 6 elements containing information about the 3D domain,
            including the number of pixels in each dimension and the conversion factors
            between pixels and real units. The first 3 elements are integers
            representing the number of pixels in the x, y, and z dimensions, respectively.
            The 4th element is a string representing the conversion factor, e.g. '10 microns/pixel'.
            The 5th element is a float representing the pixel-to-unit conversion factor, e.g. 10.0.
            The 6th element is a string representing the units, e.g. 'microns'.
        max_volume : int
            maximum allowed volume of the domain in pixels. Default is 500^3 pixels.

    Returns:
        new_dims: tuple
            typle of 6 elements containing the new dimensions of the domain,
            after decreasing its size if necessary. The first 3 elements are integers
            representing the number of pixels in the x, y, and z dimensions, respectively.
            The 4th element is a string representing the conversion factor, e.g. '10 microns/pixel'.
            The 5th element is a float representing the pixel-to-unit conversion factor, e.g. 10.0.
            The 6th element is a string representing the units, e.g. 'microns'.
        truncated : bool
            indicates whether the domain was truncated. If True, a warning message is issued.
    """
    is_2D = ccdims[6]
    default_side = round(max_volume ** (1 / 3)) if not is_2D else round(max_volume ** (1 / 2))
    old_dims = [ccdims[0], ccdims[1], ccdims[2]]
    old_volume = ccdims[0] * ccdims[1] * ccdims[2]
    if old_volume < max_volume:
        return ccdims, False
    # new_dims = ccdims
    if is_2D and old_dims[0] == old_dims[1]:
        new_dims = [default_side, default_side, 1]
    elif old_dims[0] == old_dims[1] == old_dims[2]:
        new_dims = [default_side, default_side, default_side]
    else:
        med_s = sum(old_dims) / len(old_dims)
        if not is_2D:
            proportions = [d / med_s for d in old_dims]
            new_dims = [int(default_side * p) for p in proportions]
        else:
            proportions = [d / med_s for d in old_dims[:-1]]
            new_dims = [int(default_side * p) for p in proportions]
            new_dims.append(1)
    new_dims.extend([ccdims[3], ccdims[4], ccdims[5], ccdims[6]])
    message = f"WARNING: Converted dimensions of simulation domain totaled > {round(max_volume ** (1 / 3))}**3 pixels." \
              f" \nWe have " \
              f"truncated " \
              f"the sides of the simulation.This may break the initial conditions as defined in " \
              f"PhysiCell.\nOld dimensions:{ccdims[0:3]}\nNew dimensions:{new_dims[0:3]}"
    warnings.warn(message)
    return new_dims, True


def get_diffusion_constants(d_elements):
    Ds = []
    for key, value in d_elements.items():
        if not value["use_steady_state"]:
            Ds.append(value["D"])
    return Ds


def reconvert_time_parameter(d_elements, cctime, max_D=50):
    """
    Converts the time parameter if the diffusion constants are too high.

    This function takes in diffusion parameters and a time unit and checks if the converted diffusion parameters are
    `>max_D`. If so, it issues a warning and reconverts the time unit to keep the diffusion parameters low. The new
    time unit and diffusion parameters are returned.

    The function first calls the get_diffusion_constants function to get the diffusion constants. If all diffusion
    solvers are set to use the steady state solver or the maximum diffusion constant is not greater than `max_D`, the
    original diffusion parameters and time unit are returned.

    If the maximum gamma value of the new diffusion parameters is `>=1`, the function increases tge time unit back by a
    proportion that keeps the maximum gamma value at 1 or below. Finally, the updated diffusion parameters and time
    unit are returned.

    Parameters:
    -----------
    d_elements : dict
        a dictionary of diffusion elements.
    cctime : tuple
        The previously converted time unit parameters
    max_D : float, optional
        Maximum diffusion constant allowed. Default 50

    Returns:
    -------
        d_elements : dict
            the reconverted dictionary of diffusion elements.
        new_cctime : tuple
            the reconverted tuple of cc3d time parameters
    """
    diffusion_constants = get_diffusion_constants(d_elements)

    if len(diffusion_constants) == 0 or not max(diffusion_constants) > max_D:
        return d_elements, cctime

    message = "WARNING: the converted diffusion parameters were very high, using them as is would result in a very " \
              "slow simulation. The translating software will reconvert the time unit in order to keep the diffusion" \
              " parameters low."
    warnings.warn(message)

    max_old_D = max(diffusion_constants)

    reduction_proportion = round(0.9 * max_D / max_old_D, 2)
    # reduction_proportion = float(f"{max_D / max_old_D:.3f}")

    # new_cctime = [min(int(cctime[0] / reduction_proportion), 10 ** 9),
    #               f'1 MCS = {cctime[2] * reduction_proportion} {cctime[1].split(" ")[-1]}',
    #               cctime[2] * reduction_proportion,
    #               cctime[3]]
    new_gammas = []
    for key in d_elements.keys():
        d_elements[key]["D"] *= reduction_proportion
        d_elements[key]["D_conv_factor"] *= reduction_proportion
        d_elements[key]["D_conv_factor_text"] = d_elements[key]["D_conv_factor_text"].split("=")[0] + " = " + \
                                                f'{d_elements[key]["D_conv_factor"]} ' + \
                                                d_elements[key]["D_conv_factor_text"].split(" ")[-1]

        d_elements[key]["gamma"] *= reduction_proportion
        d_elements[key]["gamma_conv_factor"] *= reduction_proportion
        d_elements[key]["gamma_conv_factor_text"] = d_elements[key]["gamma_conv_factor_text"].split("=")[0] + " = " + \
                                                    f'{d_elements[key]["gamma_conv_factor"]} ' + \
                                                    d_elements[key]["gamma_conv_factor_text"].split(" ")[-1]

        new_gammas.append(d_elements[key]["gamma"])

    new_cctime = [min(int(cctime[0] / reduction_proportion), 10 ** 9),
                  f'1 MCS = {cctime[2] * reduction_proportion} {cctime[1].split(" ")[-1]}',
                  cctime[2] * reduction_proportion,
                  cctime[3]]
    return d_elements, new_cctime

    # old_reduction_proportion = reduction_proportion
    # reduction_proportion = max(new_gammas)
    #
    # for key in d_elements.keys():
    #     d_elements[key]["D"] /= reduction_proportion
    #     d_elements[key]["D_conv_factor"] /= reduction_proportion
    #     d_elements[key]["D_conv_factor_text"] = d_elements[key]["D_conv_factor_text"].split("=")[0] + " = " + \
    #                                             f'{d_elements[key]["D_conv_factor"]} ' + \
    #                                             d_elements[key]["D_conv_factor_text"].split(" ")[-1]
    #
    #     d_elements[key]["gamma"] /= reduction_proportion
    #     d_elements[key]["gamma_conv_factor"] /= reduction_proportion
    #     d_elements[key]["gamma_conv_factor_text"] = d_elements[key]["gamma_conv_factor_text"].split("=")[0] + " = " + \
    #                                                 f'{d_elements[key]["gamma_conv_factor"]} ' + \
    #                                                 d_elements[key]["gamma_conv_factor_text"].split(" ")[-1]
    #
    # new_cctime = [min(int(cctime[0] * reduction_proportion / old_reduction_proportion), 10 ** 9),
    #               f'1 MCS = {cctime[2] * old_reduction_proportion / reduction_proportion} {cctime[1].split(" ")[-1]}',
    #               cctime[2] * reduction_proportion,
    #               cctime[3]]
    #
    # return d_elements, new_cctime
