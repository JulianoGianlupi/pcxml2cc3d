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
# has a similar implementation of phenotypes. You should install PhenoCellPy to translate the phenotypes from PhysiCell.
# Then change the default path used below with your PhenoCellPy's
# installation directory
# sys.path.extend([r'C:\modeling\phd\PhenoCellPy'])
sys.path.extend([r'D:\modeling\PhenoCellPy'])
global pcp_imp
pcp_imp = False
try:
    import PhenoCellPy as pcp

    pcp_imp = True
except BaseException:
    pass

user_data = {
    'random_seed': {
        '@type': 'int',
        '@hidden': 'true',
        '@units': 'dimensionless',
        '#text': '0'},
    'div_tum': {
        '@type': 'divider',
        '@description': '---Initial Parameters---'},
    'seeding_method': {
        '@type': 'int',
        '@units': 'micrometer',
        '@description': 'For seeding one cell => 1, For tissue seeding => 2',
        '#text': '1'},
    'div_rates': {
        '@type': 'divider',
        '@description': '---Transition Rates---'},
    'r01': {
        '@type': 'double',
        '@units': '1/min',
        '@description': 'transition rate between G0/G1 and S',
        '#text': '0.01666'},
    'r01_fixed_duration': {
        '@type': 'bool',
        '@units': '',
        '@description': 'True if transition rate (r01) are fixed duration.',
        '#text': 'false'},
    'r12': {
        '@type': 'double',
        '@units': '1/min',
        '@description': 'transition rate between S and G2',
        '#text': '0.01666'},
    'r12_fixed_duration': {
        '@type': 'bool',
        '@units': '',
        '@description': 'True if transition rate (r12) are fixed duration.',
        '#text': 'false'},
    'r23': {
        '@type': 'double',
        '@units': '1/min',
        '@description': 'transition rate between G2 and M',
        '#text': '0.01666'},
    'r23_fixed_duration': {
        '@type': 'bool',
        '@units': '',
        '@description': 'True if transition rate (r23) are fixed duration.',
        '#text': 'false'},
    'r30': {
        '@type': 'double',
        '@units': '1/min',
        '@description': 'transition rate between M and G0/G1',
        '#text': '0.01666'},
    'r30_fixed_duration': {
        '@type': 'bool',
        '@units': '',
        '@description': 'True if transition rate (r30) are fixed duration.',
        '#text': 'false'},
    'div_G0G1': {
        '@type': 'divider',
        '@description': '---Arresting phase link between G0/G1 and S---'},
    'r01_arrest': {
        '@type': 'bool',
        '@description': 'if true and conditions were matched, it arrests phase link.',
        '#text': 'false'},
    'oxygen_threshold': {
        '@type': 'double',
        '@units': 'mmHg',
        '@description': 'environmental oxygen threshold for phase transition. If less than threshold, it stops',
        '#text': '24'},
    'oxygen_gradient': {
        '@type': 'bool',
        '@units': '',
        '@description': 'True, if environment has oxygen gradient',
        '#text': 'false'},
    'div_S': {
        '@type': 'divider',
        '@description': '---Arresting phase link between S and G2---'},
    'r12_arrest': {
        '@type': 'bool',
        '@description': 'arrest function according to chemokine level',
        '#text': 'false'},
    'chemokine_threshold': {
        '@type': 'double',
        '@units': 'mM',
        '@description': 'environmental chemokine threshold for phase transition. If more, it stops',
        '#text': '10.0'},
    'div_G2': {
        '@type': 'divider',
        '@description': '---Arresting phase link between G2 an M---'},
    'r23_arrest': {
        '@type': 'bool',
        '@description': 'arrest function according to pressure level',
        '#text': 'false'},
    'pressure_threshold': {
        '@type': 'double',
        '@units': 'dimensionless',
        '@description': 'cell pressure threshold to stall this phase transition. If more, it stops',
        '#text': '1.0'},
    'div_M': {
        '@type': 'divider',
        '@description': '---Arresting phase link between M an G0/G1 (division)---'},
    'r30_arrest': {
        '@type': 'bool',
        '@description': 'arrest function according to cell volume',
        '#text': 'false'},
    'volume_threshold': {
        '@type': 'double',
        '@units': 'micron^3',
        '@description': 'target cell volume threshold. If less, it does not cycle.',
        '#text': '2490.0'}}


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
            dt = 1 / self.mcs_to_time
            self.phenotypes['CELL'] = {}
            phenotype = pcp.get_phenotype_by_name('Flow Cytometry Basic')
            self.phenotypes['CELL']['Flow Cytometry Basic'] = phenotype(dt=dt)

        cell = self.new_cell(self.CELL)
        self.cell_field[self.dim.x // 2:self.dim.x // 2 + 2, self.dim.y // 2:self.dim.y // 2 + 2, 0] = cell

        for cell in self.cell_list_by_type(self.CELL):
            cell.dict['volume'] = {'volume (None)': None, 'volume (pixels)': 32}
            cell.targetVolume = 32
            # NOTE: PC does not have an equivalent parameter, you have to
            # adjust it:
            cell.lambdaVolume = 16
            cell.dict['mechanics'] = None
            # NOTE: you are responsible for finding how this datais used in the original model
            # and re-implementing in CC3D
            cell.dict['custom_data'] = None
            if pcp_imp:
                cell.dict['phenotypes'] = self.phenotypes['CELL']
                cell.dict['current_phenotype'] = cell.dict['phenotypes']['Flow Cytometry Basic'].copy(
                )
                cell.dict['volume_conversion'] = cell.targetVolume / \
                                                 cell.dict['current_phenotype'].current_phase.volume.total
                print(cell.dict['current_phenotype'].current_phase.volume.total)
            cell.dict['phenotypes_names'] = ['Flow Cytometry Basic']

        self.shared_steppable_vars['constraints'] = self


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
            for cell in self.cell_list_by_type(self.CELL):
                # WARNING: currently you are responsible for implementing what should happen for each of
                # the flags
                changed_phase, should_be_removed, divides = \
                    cell.dict['current_phenotype'].time_step_phenotype()
                # print(changed_phase, should_be_removed, divides)
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

        # reducing parent target volume
        converted_volume = self.parent_cell.dict['volume_conversion'] * \
                           self.parent_cell.dict["current_phenotype"].current_phase.volume.total
        # print(self.parent_cell.targetVolume, converted_volume)
        self.parent_cell.targetVolume = converted_volume

        self.clone_parent_2_child()
        self.child_cell.dict["current_phenotype"] = self.parent_cell.dict["current_phenotype"].copy()

        self.parent_cell.dict["current_phenotype"].current_phase.volume.target_cytoplasm = self.parent_cell.targetVolume
        self.parent_cell.dict["current_phenotype"].current_phase.volume.cytoplasm_fluid = self.parent_cell.targetVolume
        self.parent_cell.dict["phase_index_plus_1"] = self.parent_cell.dict["current_phenotype"].current_phase.index + 1

        self.child_cell.dict["current_phenotype"].current_phase.volume.target_cytoplasm = self.parent_cell.targetVolume
        self.child_cell.dict["current_phenotype"].current_phase.volume.cytoplasm_fluid = self.parent_cell.targetVolume

        self.child_cell.dict["phase_index_plus_1"] = self.child_cell.dict["current_phenotype"].current_phase.index + 1
        self.child_cell.dict["current_phenotype"].time_in_phenotype = 0

    def finish(self):
        """
        Called after the last MCS to wrap up the simulation. Good place to close files and do post-processing
        """

    def on_stop(self):
        """
        Called if the simulation is stopped before the last MCS
        """
        self.finish()
