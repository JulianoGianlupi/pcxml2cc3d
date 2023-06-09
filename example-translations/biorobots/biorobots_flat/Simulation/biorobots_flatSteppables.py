"""
******************************************************************************************
PLEASE READ BEFORE RUNNING:
The translator program does a lot of the conversion process of a Physicell simulation into
a CompuCell3D simulation. However, just as in CC3D, the user can define several custom modules
functions, and initial conditions for their PhysiCell model. Translating custom c++ code is
beyond what an automated program can do (AI-ML language models not withstanding). On top
of that fact there are several PhysiCell concepts that do not exist in CC3D, and vice-versa.
Therefore: 1. please read *all comments* the script has placed in the converted files.
2. You are responsible for creating the initial conditions. Some are defined in csv files
others are creating with c++.
3. You are responsible for finding where <user_parameters> is used in the Physicell model,
and using it in the CC3D model
4. You are responsible for using chemical field data (e.g., chemotaxis)
5. You are responsible for defining the use of custom_data for each cell type. Find where
it is used in PhysiCell and define its use in the CC3D simulation
******************************************************************************************
"""
from cc3d.cpp.PlayerPython import *
from cc3d import CompuCellSetup
from cc3d.core.PySteppables import *
import numpy as np

import sys

# IMPORTANT: PhysiCell has a concept of cell phenotype, PhenoCellPy (https://github.com/JulianoGianlupi/PhenoCellPy)
# has a similar implementation of phenotypes. You should install PhenoCellPy to translate the Phenotypes from PhysiCell.
# Then change the default path used below with your PhenoCellPy's
# installation directory
sys.path.extend(['C:\\PhenoCellPy'])
global pcp_imp
pcp_imp = False
try:
    import PhenoCellPy as pcp
    pcp_imp = True
except BaseException:
    pass


user_data = OrderedDict([('random_seed',
                          OrderedDict([('@type',
                                        'int'),
                                       ('@units',
                                        'dimensionless'),
                                       ('#text',
                                        '0')])),
                         ('cargo_signal_D',
                          OrderedDict([('@type',
                                        'double'),
                                       ('@units',
                                        'micron/min^2'),
                                       ('#text',
                                        '1e3')])),
                         ('cargo_signal_decay',
                          OrderedDict([('@type',
                                        'double'),
                                       ('@units',
                                        '1/min'),
                                       ('#text',
                                        '.4')])),
                         ('director_signal_D',
                          OrderedDict([('@type',
                                        'double'),
                                       ('@units',
                                        'micron/min^2'),
                                       ('#text',
                                        '1e3')])),
                         ('director_signal_decay',
                          OrderedDict([('@type',
                                        'double'),
                                       ('@units',
                                        '1/min'),
                                       ('#text',
                                        '.1')])),
                         ('elastic_coefficient',
                          OrderedDict([('@type',
                                        'double'),
                                       ('@units',
                                        '1/min'),
                                       ('#text',
                                        '0.05')])),
                         ('attached_worker_migration_bias',
                          OrderedDict([('@type',
                                        'double'),
                                       ('@units',
                                        'dimensionless'),
                                       ('#text',
                                        '1.0')])),
                         ('unattached_worker_migration_bias',
                          OrderedDict([('@type',
                                        'double'),
                                       ('@units',
                                        'dimensionless'),
                                       ('#text',
                                        '0.5')])),
                         ('number_of_directors',
                          OrderedDict([('@type',
                                        'int'),
                                       ('@units',
                                        'none'),
                                       ('#text',
                                        '15')])),
                         ('number_of_cargo_clusters',
                          OrderedDict([('@type',
                                        'int'),
                                       ('@units',
                                        'none'),
                                       ('#text',
                                        '100')])),
                         ('number_of_workers',
                          OrderedDict([('@type',
                                        'int'),
                                       ('@units',
                                        'none'),
                                       ('#text',
                                        '50')])),
                         ('drop_threshold',
                          OrderedDict([('@type',
                                        'double'),
                                       ('@units',
                                        'dimensionless'),
                                       ('#text',
                                        '0.4')])),
                         ('worker_color',
                          OrderedDict([('@type',
                                        'string'),
                                       ('@units',
                                        'none'),
                                       ('#text',
                                        'red')])),
                         ('cargo_color',
                          OrderedDict([('@type',
                                        'string'),
                                       ('@units',
                                        'none'),
                                       ('#text',
                                        'blue')])),
                         ('director_color',
                          OrderedDict([('@type',
                                        'string'),
                                       ('@units',
                                        'none'),
                                       ('#text',
                                        'limegreen')]))])


class ConstraintsSteppable(SteppableBasePy):

    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        """
        Called before MCS=0 while building the initial simulation
        """
        self.pixel_to_space = float(self.get_xml_element(
            'pixel_to_space').cdata)  # pixel/[unit], see xml for units
        self.mcs_to_time = float(self.get_xml_element(
            'mcs_to_time').cdata)  # MCS/[unit], see xml for units

        if pcp_imp:
            self.phenotypes = {}
            dt = 1/self.mcs_to_time
            self.phenotypes['default'] = {}
            phenotype = pcp.get_phenotype_by_name('Simple Live')
            self.phenotypes['default']['Simple Live'] = phenotype(dt=dt,
                                                                  time_unit='min',
                                                                  fixed_durations=[
                                                                      False],
                                                                  phase_durations=[
                                                                      9e+99],
                                                                  cytoplasm_volume_change_rate=[
                                                                      0.0045],
                                                                  nuclear_volume_change_rate=[
                                                                      0.0055],
                                                                  calcification_rate=[
                                                                      0.0],
                                                                  calcified_fraction=[
                                                                      0.0],
                                                                  target_fluid_fraction=[
                                                                      0.75],
                                                                  nuclear_fluid=[
                                                                      405.0],
                                                                  nuclear_solid=[
                                                                      135.0],
                                                                  nuclear_solid_target=[
                                                                      135.0],
                                                                  cytoplasm_fluid=[
                                                                      1465.5],
                                                                  cytoplasm_solid=[
                                                                      488.5],
                                                                  cytoplasm_solid_target=[
                                                                      488.5],
                                                                  target_cytoplasm_to_nuclear_ratio=[
                                                                      3.6185185185185187],
                                                                  fluid_change_rate=[0.05])
            phenotype = pcp.get_phenotype_by_name('Standard apoptosis model')
            self.phenotypes['default']['Standard apoptosis model'] = phenotype(dt=dt,
                                                                               time_unit='min',
                                                                               fixed_durations=[
                                                                                   True],
                                                                               phase_durations=[
                                                                                   516.0],
                                                                               cytoplasm_volume_change_rate=[
                                                                                   0.0166667],
                                                                               nuclear_volume_change_rate=[
                                                                                   0.00583333],
                                                                               calcification_rate=[
                                                                                   0.0],
                                                                               calcified_fraction=[
                                                                                   0],
                                                                               target_fluid_fraction=[
                                                                                   0.75],
                                                                               nuclear_fluid=[
                                                                                   None],
                                                                               nuclear_solid=[
                                                                                   None],
                                                                               nuclear_solid_target=[
                                                                                   None],
                                                                               cytoplasm_fluid=[
                                                                                   None],
                                                                               cytoplasm_solid=[
                                                                                   None],
                                                                               cytoplasm_solid_target=[
                                                                                   None],
                                                                               target_cytoplasm_to_nuclear_ratio=[
                                                                                   None],
                                                                               fluid_change_rate=[0.05])
            phenotype = pcp.get_phenotype_by_name('Standard necrosis model')
            self.phenotypes['default']['Standard necrosis model'] = phenotype(dt=dt,
                                                                              time_unit='min',
                                                                              fixed_durations=[
                                                                                  True, True],
                                                                              phase_durations=[
                                                                                  9e+99, 86400.0],
                                                                              cytoplasm_volume_change_rate=[
                                                                                  0.0166667, 0.0166667],
                                                                              nuclear_volume_change_rate=[
                                                                                  0.00583333, 0.00583333],
                                                                              calcification_rate=[
                                                                                  0.0, 0.0],
                                                                              calcified_fraction=[
                                                                                  0, 0],
                                                                              target_fluid_fraction=[
                                                                                  0.75, 0.75],
                                                                              nuclear_fluid=[
                                                                                  None, None],
                                                                              nuclear_solid=[
                                                                                  None, None],
                                                                              nuclear_solid_target=[
                                                                                  None, None],
                                                                              cytoplasm_fluid=[
                                                                                  None, None],
                                                                              cytoplasm_solid=[
                                                                                  None, None],
                                                                              cytoplasm_solid_target=[
                                                                                  None, None],
                                                                              target_cytoplasm_to_nuclear_ratio=[
                                                                                  None, None],
                                                                              fluid_change_rate=[0.05, 0.0])
            dt = 1/self.mcs_to_time
            self.phenotypes['director_cell'] = {}
            phenotype = pcp.get_phenotype_by_name('Simple Live')
            self.phenotypes['director_cell']['Simple Live'] = phenotype(dt=dt,
                                                                        time_unit='min',
                                                                        fixed_durations=[
                                                                            False],
                                                                        phase_durations=[
                                                                            9e+99],
                                                                        cytoplasm_volume_change_rate=[
                                                                            0.0045],
                                                                        nuclear_volume_change_rate=[
                                                                            0.0055],
                                                                        calcification_rate=[
                                                                            0.0],
                                                                        calcified_fraction=[
                                                                            0.0],
                                                                        target_fluid_fraction=[
                                                                            0.75],
                                                                        nuclear_fluid=[
                                                                            405.0],
                                                                        nuclear_solid=[
                                                                            135.0],
                                                                        nuclear_solid_target=[
                                                                            135.0],
                                                                        cytoplasm_fluid=[
                                                                            1465.5],
                                                                        cytoplasm_solid=[
                                                                            488.5],
                                                                        cytoplasm_solid_target=[
                                                                            488.5],
                                                                        target_cytoplasm_to_nuclear_ratio=[
                                                                            3.6185185185185187],
                                                                        fluid_change_rate=[0.05])
            phenotype = pcp.get_phenotype_by_name('Standard apoptosis model')
            self.phenotypes['director_cell']['Standard apoptosis model'] = phenotype(dt=dt,
                                                                                     time_unit='min',
                                                                                     fixed_durations=[
                                                                                         True],
                                                                                     phase_durations=[
                                                                                         516.0],
                                                                                     cytoplasm_volume_change_rate=[
                                                                                         0.0166667],
                                                                                     nuclear_volume_change_rate=[
                                                                                         0.00583333],
                                                                                     calcification_rate=[
                                                                                         0.0],
                                                                                     calcified_fraction=[
                                                                                         0],
                                                                                     target_fluid_fraction=[
                                                                                         0.75],
                                                                                     nuclear_fluid=[
                                                                                         None],
                                                                                     nuclear_solid=[
                                                                                         None],
                                                                                     nuclear_solid_target=[
                                                                                         None],
                                                                                     cytoplasm_fluid=[
                                                                                         None],
                                                                                     cytoplasm_solid=[
                                                                                         None],
                                                                                     cytoplasm_solid_target=[
                                                                                         None],
                                                                                     target_cytoplasm_to_nuclear_ratio=[
                                                                                         None],
                                                                                     fluid_change_rate=[0.05])
            phenotype = pcp.get_phenotype_by_name('Standard necrosis model')
            self.phenotypes['director_cell']['Standard necrosis model'] = phenotype(dt=dt,
                                                                                    time_unit='min',
                                                                                    fixed_durations=[
                                                                                        True, True],
                                                                                    phase_durations=[
                                                                                        9e+99, 86400.0],
                                                                                    cytoplasm_volume_change_rate=[
                                                                                        0.0166667, 0.0166667],
                                                                                    nuclear_volume_change_rate=[
                                                                                        0.00583333, 0.00583333],
                                                                                    calcification_rate=[
                                                                                        0.0, 0.0],
                                                                                    calcified_fraction=[
                                                                                        0, 0],
                                                                                    target_fluid_fraction=[
                                                                                        0.75, 0.75],
                                                                                    nuclear_fluid=[
                                                                                        None, None],
                                                                                    nuclear_solid=[
                                                                                        None, None],
                                                                                    nuclear_solid_target=[
                                                                                        None, None],
                                                                                    cytoplasm_fluid=[
                                                                                        None, None],
                                                                                    cytoplasm_solid=[
                                                                                        None, None],
                                                                                    cytoplasm_solid_target=[
                                                                                        None, None],
                                                                                    target_cytoplasm_to_nuclear_ratio=[
                                                                                        None, None],
                                                                                    fluid_change_rate=[0.05, 0.0])
            dt = 1/self.mcs_to_time
            self.phenotypes['cargo_cell'] = {}
            phenotype = pcp.get_phenotype_by_name('Simple Live')
            self.phenotypes['cargo_cell']['Simple Live'] = phenotype(dt=dt,
                                                                     time_unit='min',
                                                                     fixed_durations=[
                                                                         False],
                                                                     phase_durations=[
                                                                         9e+99],
                                                                     cytoplasm_volume_change_rate=[
                                                                         0.0045],
                                                                     nuclear_volume_change_rate=[
                                                                         0.0055],
                                                                     calcification_rate=[
                                                                         0.0],
                                                                     calcified_fraction=[
                                                                         0.0],
                                                                     target_fluid_fraction=[
                                                                         0.75],
                                                                     nuclear_fluid=[
                                                                         405.0],
                                                                     nuclear_solid=[
                                                                         135.0],
                                                                     nuclear_solid_target=[
                                                                         135.0],
                                                                     cytoplasm_fluid=[
                                                                         1465.5],
                                                                     cytoplasm_solid=[
                                                                         488.5],
                                                                     cytoplasm_solid_target=[
                                                                         488.5],
                                                                     target_cytoplasm_to_nuclear_ratio=[
                                                                         3.6185185185185187],
                                                                     fluid_change_rate=[0.05])
            phenotype = pcp.get_phenotype_by_name('Standard apoptosis model')
            self.phenotypes['cargo_cell']['Standard apoptosis model'] = phenotype(dt=dt,
                                                                                  time_unit='min',
                                                                                  fixed_durations=[
                                                                                      True],
                                                                                  phase_durations=[
                                                                                      516.0],
                                                                                  cytoplasm_volume_change_rate=[
                                                                                      0.0166667],
                                                                                  nuclear_volume_change_rate=[
                                                                                      0.00583333],
                                                                                  calcification_rate=[
                                                                                      0.0],
                                                                                  calcified_fraction=[
                                                                                      0],
                                                                                  target_fluid_fraction=[
                                                                                      0.75],
                                                                                  nuclear_fluid=[
                                                                                      None],
                                                                                  nuclear_solid=[
                                                                                      None],
                                                                                  nuclear_solid_target=[
                                                                                      None],
                                                                                  cytoplasm_fluid=[
                                                                                      None],
                                                                                  cytoplasm_solid=[
                                                                                      None],
                                                                                  cytoplasm_solid_target=[
                                                                                      None],
                                                                                  target_cytoplasm_to_nuclear_ratio=[
                                                                                      None],
                                                                                  fluid_change_rate=[0.05])
            phenotype = pcp.get_phenotype_by_name('Standard necrosis model')
            self.phenotypes['cargo_cell']['Standard necrosis model'] = phenotype(dt=dt,
                                                                                 time_unit='min',
                                                                                 fixed_durations=[
                                                                                     True, True],
                                                                                 phase_durations=[
                                                                                     9e+99, 86400.0],
                                                                                 cytoplasm_volume_change_rate=[
                                                                                     0.0166667, 0.0166667],
                                                                                 nuclear_volume_change_rate=[
                                                                                     0.00583333, 0.00583333],
                                                                                 calcification_rate=[
                                                                                     0.0, 0.0],
                                                                                 calcified_fraction=[
                                                                                     0, 0],
                                                                                 target_fluid_fraction=[
                                                                                     0.75, 0.75],
                                                                                 nuclear_fluid=[
                                                                                     None, None],
                                                                                 nuclear_solid=[
                                                                                     None, None],
                                                                                 nuclear_solid_target=[
                                                                                     None, None],
                                                                                 cytoplasm_fluid=[
                                                                                     None, None],
                                                                                 cytoplasm_solid=[
                                                                                     None, None],
                                                                                 cytoplasm_solid_target=[
                                                                                     None, None],
                                                                                 target_cytoplasm_to_nuclear_ratio=[
                                                                                     None, None],
                                                                                 fluid_change_rate=[0.05, 0.0])
            dt = 1/self.mcs_to_time
            self.phenotypes['worker_cell'] = {}
            phenotype = pcp.get_phenotype_by_name('Simple Live')
            self.phenotypes['worker_cell']['Simple Live'] = phenotype(dt=dt,
                                                                      time_unit='min',
                                                                      fixed_durations=[
                                                                          False],
                                                                      phase_durations=[
                                                                          9e+99],
                                                                      cytoplasm_volume_change_rate=[
                                                                          0.0045],
                                                                      nuclear_volume_change_rate=[
                                                                          0.0055],
                                                                      calcification_rate=[
                                                                          0.0],
                                                                      calcified_fraction=[
                                                                          0.0],
                                                                      target_fluid_fraction=[
                                                                          0.75],
                                                                      nuclear_fluid=[
                                                                          405.0],
                                                                      nuclear_solid=[
                                                                          135.0],
                                                                      nuclear_solid_target=[
                                                                          135.0],
                                                                      cytoplasm_fluid=[
                                                                          1465.5],
                                                                      cytoplasm_solid=[
                                                                          488.5],
                                                                      cytoplasm_solid_target=[
                                                                          488.5],
                                                                      target_cytoplasm_to_nuclear_ratio=[
                                                                          3.6185185185185187],
                                                                      fluid_change_rate=[0.05])
            phenotype = pcp.get_phenotype_by_name('Standard apoptosis model')
            self.phenotypes['worker_cell']['Standard apoptosis model'] = phenotype(dt=dt,
                                                                                   time_unit='min',
                                                                                   fixed_durations=[
                                                                                       True],
                                                                                   phase_durations=[
                                                                                       516.0],
                                                                                   cytoplasm_volume_change_rate=[
                                                                                       0.0166667],
                                                                                   nuclear_volume_change_rate=[
                                                                                       0.00583333],
                                                                                   calcification_rate=[
                                                                                       0.0],
                                                                                   calcified_fraction=[
                                                                                       0],
                                                                                   target_fluid_fraction=[
                                                                                       0.75],
                                                                                   nuclear_fluid=[
                                                                                       None],
                                                                                   nuclear_solid=[
                                                                                       None],
                                                                                   nuclear_solid_target=[
                                                                                       None],
                                                                                   cytoplasm_fluid=[
                                                                                       None],
                                                                                   cytoplasm_solid=[
                                                                                       None],
                                                                                   cytoplasm_solid_target=[
                                                                                       None],
                                                                                   target_cytoplasm_to_nuclear_ratio=[
                                                                                       None],
                                                                                   fluid_change_rate=[0.05])
            phenotype = pcp.get_phenotype_by_name('Standard necrosis model')
            self.phenotypes['worker_cell']['Standard necrosis model'] = phenotype(dt=dt,
                                                                                  time_unit='min',
                                                                                  fixed_durations=[
                                                                                      True, True],
                                                                                  phase_durations=[
                                                                                      9e+99, 86400.0],
                                                                                  cytoplasm_volume_change_rate=[
                                                                                      0.0166667, 0.0166667],
                                                                                  nuclear_volume_change_rate=[
                                                                                      0.00583333, 0.00583333],
                                                                                  calcification_rate=[
                                                                                      0.0, 0.0],
                                                                                  calcified_fraction=[
                                                                                      0, 0],
                                                                                  target_fluid_fraction=[
                                                                                      0.75, 0.75],
                                                                                  nuclear_fluid=[
                                                                                      None, None],
                                                                                  nuclear_solid=[
                                                                                      None, None],
                                                                                  nuclear_solid_target=[
                                                                                      None, None],
                                                                                  cytoplasm_fluid=[
                                                                                      None, None],
                                                                                  cytoplasm_solid=[
                                                                                      None, None],
                                                                                  cytoplasm_solid_target=[
                                                                                      None, None],
                                                                                  target_cytoplasm_to_nuclear_ratio=[
                                                                                      None, None],
                                                                                  fluid_change_rate=[0.05, 0.0])

        for cell in self.cell_list_by_type(self.DEFAULT):
            cell.dict['volume'] = {
                'volume (micron^3)': 2494.0,
                'volume (pixels)': 8.105500000000003}
            cell.targetVolume = 8.105500000000003
            # NOTE: PC does not have an equivalent parameter, you have to
            # adjust it:
            cell.lambdaVolume = 8
            cell.dict['mechanics'] = {
                'cell_cell_adhesion_strength': {
                    'units': 'micron/min',
                    'value': 0.4},
                'cell_cell_repulsion_strength': {
                    'units': 'micron/min',
                    'value': 10.0},
                'relative_maximum_adhesion_distance': {
                    'units': 'dimensionless',
                    'value': 1.25}}
            # NOTE: you are responsible for finding how this datais used in the original model
            # and re-implementing in CC3D
            cell.dict['custom_data'] = OrderedDict(
                [('receptor', OrderedDict([('@units', 'dimensionless'), ('#text', '0.0')]))])
            if pcp_imp:
                cell.dict['phenotypes'] = self.phenotypes['default']
                cell.dict['current_phenotype'] = cell.dict['phenotypes']['Simple Live'].copy(
                )
                cell.dict['volume_conversion'] = cell.targetVolume / \
                    cell.dict['current_phenotype'].current_phase.volume.total
            cell.dict['phenotypes_names'] = [
                'Simple Live',
                'Standard apoptosis model',
                'Standard necrosis model']
            cell.dict['director_signal'] = {
                'secretion_rate': 0.0,
                'secretion_unit': '1/min',
                'secretion_target': 1.0,
                'uptake_rate': 0.0,
                'uptake_unit': '1/min',
                'net_export': 0.0,
                'net_export_unit': 'total substrate/min',
                'secretion_rate_MCS': 0.0,
                'net_export_MCS': 0.0,
                'uptake_rate_MCS': 0.0}
            cell.dict['cargo_signal'] = {
                'secretion_rate': 0.0,
                'secretion_unit': '1/min',
                'secretion_target': 1.0,
                'uptake_rate': 0.0,
                'uptake_unit': '1/min',
                'net_export': 0.0,
                'net_export_unit': 'total substrate/min',
                'secretion_rate_MCS': 0.0,
                'net_export_MCS': 0.0,
                'uptake_rate_MCS': 0.0}

        for cell in self.cell_list_by_type(self.DIRECTOR_CELL):
            cell.dict['volume'] = {
                'volume (micron^3)': 2494.0,
                'volume (pixels)': 8.105500000000003}
            cell.targetVolume = 8.105500000000003
            # NOTE: PC does not have an equivalent parameter, you have to
            # adjust it:
            cell.lambdaVolume = 8
            cell.dict['mechanics'] = {
                'cell_cell_adhesion_strength': {
                    'units': 'micron/min',
                    'value': 0.4},
                'cell_cell_repulsion_strength': {
                    'units': 'micron/min',
                    'value': 10.0},
                'relative_maximum_adhesion_distance': {
                    'units': 'dimensionless',
                    'value': 1.25}}
            # NOTE: you are responsible for finding how this datais used in the original model
            # and re-implementing in CC3D
            cell.dict['custom_data'] = OrderedDict(
                [('receptor', OrderedDict([('@units', 'dimensionless'), ('#text', '0.0')]))])
            if pcp_imp:
                cell.dict['phenotypes'] = self.phenotypes['director_cell']
                cell.dict['current_phenotype'] = cell.dict['phenotypes']['Simple Live'].copy(
                )
                cell.dict['volume_conversion'] = cell.targetVolume / \
                    cell.dict['current_phenotype'].current_phase.volume.total
            cell.dict['phenotypes_names'] = [
                'Simple Live',
                'Standard apoptosis model',
                'Standard necrosis model']
            cell.dict['director_signal'] = {
                'secretion_rate': 9.9,
                'secretion_unit': '1/min',
                'secretion_target': 1.0,
                'uptake_rate': 0.0,
                'uptake_unit': '1/min',
                'net_export': 0.0,
                'net_export_unit': 'total substrate/min',
                'secretion_rate_MCS': 3.6666666666666665,
                'net_export_MCS': 0.0,
                'uptake_rate_MCS': 0.0}
            cell.dict['cargo_signal'] = {
                'secretion_rate': 0.0,
                'secretion_unit': '1/min',
                'secretion_target': 1.0,
                'uptake_rate': 0.0,
                'uptake_unit': '1/min',
                'net_export': 0.0,
                'net_export_unit': 'total substrate/min',
                'secretion_rate_MCS': 0.0,
                'net_export_MCS': 0.0,
                'uptake_rate_MCS': 0.0}

        for cell in self.cell_list_by_type(self.CARGO_CELL):
            cell.dict['volume'] = {
                'volume (micron^3)': 2494.0,
                'volume (pixels)': 8.105500000000003}
            cell.targetVolume = 8.105500000000003
            # NOTE: PC does not have an equivalent parameter, you have to
            # adjust it:
            cell.lambdaVolume = 8
            cell.dict['mechanics'] = {
                'cell_cell_adhesion_strength': {
                    'units': 'micron/min',
                    'value': 0.4},
                'cell_cell_repulsion_strength': {
                    'units': 'micron/min',
                    'value': 10.0},
                'relative_maximum_adhesion_distance': {
                    'units': 'dimensionless',
                    'value': 1.25}}
            # NOTE: you are responsible for finding how this datais used in the original model
            # and re-implementing in CC3D
            cell.dict['custom_data'] = OrderedDict(
                [('receptor', OrderedDict([('@units', 'dimensionless'), ('#text', '1.0')]))])
            if pcp_imp:
                cell.dict['phenotypes'] = self.phenotypes['cargo_cell']
                cell.dict['current_phenotype'] = cell.dict['phenotypes']['Simple Live'].copy(
                )
                cell.dict['volume_conversion'] = cell.targetVolume / \
                    cell.dict['current_phenotype'].current_phase.volume.total
            cell.dict['phenotypes_names'] = [
                'Simple Live',
                'Standard apoptosis model',
                'Standard necrosis model']
            cell.dict['director_signal'] = {
                'secretion_rate': 0.0,
                'secretion_unit': '1/min',
                'secretion_target': 1.0,
                'uptake_rate': 0.0,
                'uptake_unit': '1/min',
                'net_export': 0.0,
                'net_export_unit': 'total substrate/min',
                'secretion_rate_MCS': 0.0,
                'net_export_MCS': 0.0,
                'uptake_rate_MCS': 0.0}
            cell.dict['cargo_signal'] = {
                'secretion_rate': 9.9,
                'secretion_unit': '1/min',
                'secretion_target': 1.0,
                'uptake_rate': 0.0,
                'uptake_unit': '1/min',
                'net_export': 0.0,
                'net_export_unit': 'total substrate/min',
                'secretion_rate_MCS': 3.6666666666666665,
                'net_export_MCS': 0.0,
                'uptake_rate_MCS': 0.0}

        for cell in self.cell_list_by_type(self.WORKER_CELL):
            cell.dict['volume'] = {
                'volume (micron^3)': 2494.0,
                'volume (pixels)': 8.105500000000003}
            cell.targetVolume = 8.105500000000003
            # NOTE: PC does not have an equivalent parameter, you have to
            # adjust it:
            cell.lambdaVolume = 8
            cell.dict['mechanics'] = {
                'cell_cell_adhesion_strength': {
                    'units': 'micron/min',
                    'value': 0.4},
                'cell_cell_repulsion_strength': {
                    'units': 'micron/min',
                    'value': 10.0},
                'relative_maximum_adhesion_distance': {
                    'units': 'dimensionless',
                    'value': 1.25}}
            # NOTE: you are responsible for finding how this datais used in the original model
            # and re-implementing in CC3D
            cell.dict['custom_data'] = OrderedDict(
                [('receptor', OrderedDict([('@units', 'dimensionless'), ('#text', '1.0')]))])
            if pcp_imp:
                cell.dict['phenotypes'] = self.phenotypes['worker_cell']
                cell.dict['current_phenotype'] = cell.dict['phenotypes']['Simple Live'].copy(
                )
                cell.dict['volume_conversion'] = cell.targetVolume / \
                    cell.dict['current_phenotype'].current_phase.volume.total
            cell.dict['phenotypes_names'] = [
                'Simple Live',
                'Standard apoptosis model',
                'Standard necrosis model']
            cell.dict['director_signal'] = {
                'secretion_rate': 0.0,
                'secretion_unit': '1/min',
                'secretion_target': 1.0,
                'uptake_rate': 0.0,
                'uptake_unit': '1/min',
                'net_export': 0.0,
                'net_export_unit': 'total substrate/min',
                'secretion_rate_MCS': 0.0,
                'net_export_MCS': 0.0,
                'uptake_rate_MCS': 0.0}
            cell.dict['cargo_signal'] = {
                'secretion_rate': 0.0,
                'secretion_unit': '1/min',
                'secretion_target': 1.0,
                'uptake_rate': 0.0,
                'uptake_unit': '1/min',
                'net_export': 0.0,
                'net_export_unit': 'total substrate/min',
                'secretion_rate_MCS': 0.0,
                'net_export_MCS': 0.0,
                'uptake_rate_MCS': 0.0}

        self.build_wall(self.WALL)
        self.shared_steppable_vars['constraints'] = self


class SecretionUptakeSteppable(SteppableBasePy):

    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        """
        Called before MCS=0 while building the initial simulation
        """
        self.pixel_to_space = float(self.get_xml_element(
            'pixel_to_space').cdata)  # pixel/[unit], see xml for units
        self.mcs_to_time = float(self.get_xml_element(
            'mcs_to_time').cdata)  # MCS/[unit], see xml for units
        self.secretors = {
            'director_signal': self.get_field_secretor('director_signal'),
            'cargo_signal': self.get_field_secretor('cargo_signal')}

    def step(self, mcs):
        """
        Called every frequency MCS while executing the simulation

        :param mcs: current Monte Carlo step
        """

        for field_name, secretor in self.secretors.items():
            for cell in self.cell_list_by_type(self.DEFAULT):
                if field_name in cell.dict.keys():
                    data = cell.dict[field_name]
                    seen = secretor.amountSeenByCell(cell)
# WARNING: PhysiCell has a concept of "target secretion" that CompuCell3D does not.
# The translating program attempts to implement it, but it may not be a 1
# to 1 conversion.
                    net_secretion = max(
                        0, data['secretion_rate_MCS'] * (data['secretion_target'] - seen)) + data['net_export_MCS']
                    # In PhysiCell cells are point-like, in CC3D they have an arbitrary shape. With this
                    # CC3D allows several different secretion locations: over the whole cell (what the translator uses),
                    # just inside the cell surface, just outside the surface,
                    # at the surface. You should explore the options
                    if net_secretion:
                        secretor.secreteInsideCell(cell, net_secretion)
                    if data['uptake_rate']:
                        secretor.uptakeInsideCell(
                            cell, 1e10, data['uptake_rate'])
            for cell in self.cell_list_by_type(self.DIRECTOR_CELL):
                if field_name in cell.dict.keys():
                    data = cell.dict[field_name]
                    seen = secretor.amountSeenByCell(cell)
# WARNING: PhysiCell has a concept of "target secretion" that CompuCell3D does not.
# The translating program attempts to implement it, but it may not be a 1
# to 1 conversion.
                    net_secretion = max(
                        0, data['secretion_rate_MCS'] * (data['secretion_target'] - seen)) + data['net_export_MCS']
                    # In PhysiCell cells are point-like, in CC3D they have an arbitrary shape. With this
                    # CC3D allows several different secretion locations: over the whole cell (what the translator uses),
                    # just inside the cell surface, just outside the surface,
                    # at the surface. You should explore the options
                    if net_secretion:
                        secretor.secreteInsideCell(cell, net_secretion)
                    if data['uptake_rate']:
                        secretor.uptakeInsideCell(
                            cell, 1e10, data['uptake_rate'])
            for cell in self.cell_list_by_type(self.CARGO_CELL):
                if field_name in cell.dict.keys():
                    data = cell.dict[field_name]
                    seen = secretor.amountSeenByCell(cell)
# WARNING: PhysiCell has a concept of "target secretion" that CompuCell3D does not.
# The translating program attempts to implement it, but it may not be a 1
# to 1 conversion.
                    net_secretion = max(
                        0, data['secretion_rate_MCS'] * (data['secretion_target'] - seen)) + data['net_export_MCS']
                    # In PhysiCell cells are point-like, in CC3D they have an arbitrary shape. With this
                    # CC3D allows several different secretion locations: over the whole cell (what the translator uses),
                    # just inside the cell surface, just outside the surface,
                    # at the surface. You should explore the options
                    if net_secretion:
                        secretor.secreteInsideCell(cell, net_secretion)
                    if data['uptake_rate']:
                        secretor.uptakeInsideCell(
                            cell, 1e10, data['uptake_rate'])
            for cell in self.cell_list_by_type(self.WORKER_CELL):
                if field_name in cell.dict.keys():
                    data = cell.dict[field_name]
                    seen = secretor.amountSeenByCell(cell)
# WARNING: PhysiCell has a concept of "target secretion" that CompuCell3D does not.
# The translating program attempts to implement it, but it may not be a 1
# to 1 conversion.
                    net_secretion = max(
                        0, data['secretion_rate_MCS'] * (data['secretion_target'] - seen)) + data['net_export_MCS']
                    # In PhysiCell cells are point-like, in CC3D they have an arbitrary shape. With this
                    # CC3D allows several different secretion locations: over the whole cell (what the translator uses),
                    # just inside the cell surface, just outside the surface,
                    # at the surface. You should explore the options
                    if net_secretion:
                        secretor.secreteInsideCell(cell, net_secretion)
                    if data['uptake_rate']:
                        secretor.uptakeInsideCell(
                            cell, 1e10, data['uptake_rate'])

    def finish(self):
        """
        Called after the last MCS to wrap up the simulation. Good place to close files and do post-processing
        """

    def on_stop(self):
        """
        Called if the simulation is stopped before the last MCS
        """
        self.finish()


class PhenotypeSteppable(MitosisSteppableBase):

    def __init__(self, frequency=1):
        MitosisSteppableBase.__init__(self, frequency)

    def start(self):
        """
        Called before MCS=0 while building the initial simulation
        """
        self.pixel_to_space = float(self.get_xml_element(
            'pixel_to_space').cdata)  # pixel/[unit], see xml for units
        self.mcs_to_time = float(self.get_xml_element(
            'mcs_to_time').cdata)  # MCS/[unit], see xml for units
        pass

    def step(self, mcs):
        """
        Called every frequency MCS while executing the simulation

        :param mcs: current Monte Carlo step
        """

        cells_to_divide = []
        if pcp_imp:
            pass
            for cell in self.cell_list_by_type(self.DEFAULT):
                # WARNING: currently you are responsible for implementing what should happen for each of
                # the flags
                changed_phase, should_be_removed, divides = \
                    cell.dict['current_phenotype'].time_step_phenotype()
                if divides:
                    cells_to_divide.append(cell)
                cell.targetVolume = cell.dict['volume_conversion'] * \
                    cell.dict['current_phenotype'].current_phase.volume.total
            for cell in self.cell_list_by_type(self.DIRECTOR_CELL):
                # WARNING: currently you are responsible for implementing what should happen for each of
                # the flags
                changed_phase, should_be_removed, divides = \
                    cell.dict['current_phenotype'].time_step_phenotype()
                if divides:
                    cells_to_divide.append(cell)
                cell.targetVolume = cell.dict['volume_conversion'] * \
                    cell.dict['current_phenotype'].current_phase.volume.total
            for cell in self.cell_list_by_type(self.CARGO_CELL):
                # WARNING: currently you are responsible for implementing what should happen for each of
                # the flags
                changed_phase, should_be_removed, divides = \
                    cell.dict['current_phenotype'].time_step_phenotype()
                if divides:
                    cells_to_divide.append(cell)
                cell.targetVolume = cell.dict['volume_conversion'] * \
                    cell.dict['current_phenotype'].current_phase.volume.total
            for cell in self.cell_list_by_type(self.WORKER_CELL):
                # WARNING: currently you are responsible for implementing what should happen for each of
                # the flags
                changed_phase, should_be_removed, divides = \
                    cell.dict['current_phenotype'].time_step_phenotype()
                if divides:
                    cells_to_divide.append(cell)
                cell.targetVolume = cell.dict['volume_conversion'] * \
                    cell.dict['current_phenotype'].current_phase.volume.total
            for cell in cells_to_divide:
                # WARNING: As cells in CC3D have shape, they can be divided
                # along their minor/major axis, randomly in half, or along a
                # specific vector
                self.divide_cell_random_orientation(cell)
                # self.divide_cell_orientation_vector_based(cell, 1, 1, 0)
                # self.divide_cell_along_major_axis(cell)
                # self.divide_cell_along_minor_axis(cell)

    def update_attributes(self):
        self.parent_cell.targetVolume /= 2.0
        self.clone_parent_2_child()

    def finish(self):
        """
        Called after the last MCS to wrap up the simulation. Good place to close files and do post-processing
        """

    def on_stop(self):
        """
        Called if the simulation is stopped before the last MCS
        """
        self.finish()
