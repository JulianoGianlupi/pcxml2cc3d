
def get_steppables_names(step_string):
    parts = step_string.split("class")

    step_names = []

    for part in parts:
        if "Steppable" in part and "import" not in part:
            name = part.split("\n")[0].split("(")[0].replace(" ", "")
            step_names.append(name)
    return step_names


if __name__ == "__main__":
    s = '''from cc3d.cpp.PlayerPython import *
from cc3d import CompuCellSetup
from cc3d.core.PySteppables import *
import numpy as np

class ConstraintsSteppable(SteppableBasePy):

	def __init__(self, frequency=1):
		SteppableBasePy.__init__(self,frequency)

	def start(self):
		"""
		Called before MCS=0 while building the initial simulation
		"""

		for cell in self.cell_list_by_type(self.CANCER_CELL):
			cell.dict['volume']={'volume (micron^3)': 2494.0, 'volume (pixels)': 0.3117500000000001}
			cell.targetVolume = 0.3117500000000001
			cell.lambdaVolume = 8 # NOTE: PC does not have an equivalent parameter. You have to adjust it
			cell.dict['mechanics']={'cell_cell_adhesion_strength': {'units': 'micron/min', 'value': 0.0}, 'cell_cell_repulsion_strength': {'units': 'micron/min', 'value': 10.0}, 'relative_maximum_adhesion_distance': {'units': 'dimensionless', 'value': 1.25}}
			cell.dict['oxygen']={'secretion_rate': 0.0, 'secretion_unit': '1/min', 'secretion_target': 38.0, 'uptake_rate': 10.0, 'uptake_unit': '1/min', 'net_export': 0.0, 'net_export_unit': 'total substrate/min', 'secretion_rate_MCS': 0.0, 'secretion_comment': '#WARNING: PhysiCell has a concept of "target secretion" that CompuCell3D does not. \n#The translating program attempts to implement it, but it may not be a 1 to 1 conversion.', 'net_export_MCS': 0.0, 'net_secretion_comment': '', 'uptake_rate_MCS': 1.0, 'uptake_comment': '#WARNING: To avoid negative concentrations, in CompuCell3D uptake is "bounded." \n# If the amount that would be uptaken is larger than the value at that pixel,\n# the uptake will be a set ratio of the amount available.\n# The conversion program uses 1 as the ratio,\n# you may want to revisit this.'}
			cell.dict['immunostimulatory_factor']={'secretion_rate': 0.0, 'secretion_unit': '1/min', 'secretion_target': 1.0, 'uptake_rate': 0.0, 'uptake_unit': '1/min', 'net_export': 0.0, 'net_export_unit': 'total substrate/min', 'secretion_rate_MCS': 0.0, 'secretion_comment': '#WARNING: PhysiCell has a concept of "target secretion" that CompuCell3D does not. \n#The translating program attempts to implement it, but it may not be a 1 to 1 conversion.', 'net_export_MCS': 0.0, 'net_secretion_comment': '', 'uptake_rate_MCS': 0.0, 'uptake_comment': '#WARNING: To avoid negative concentrations, in CompuCell3D uptake is "bounded." \n# If the amount that would be uptaken is larger than the value at that pixel,\n# the uptake will be a set ratio of the amount available.\n# The conversion program uses 1 as the ratio,\n# you may want to revisit this.'}


		for cell in self.cell_list_by_type(self.IMMUNE_CELL):
			cell.dict['volume']={'volume (micron^3)': 2494.0, 'volume (pixels)': 0.3117500000000001}
			cell.targetVolume = 0.3117500000000001
			cell.lambdaVolume = 8 # NOTE: PC does not have an equivalent parameter. You have to adjust it
			cell.dict['mechanics']={'cell_cell_adhesion_strength': {'units': 'micron/min', 'value': 0.0}, 'cell_cell_repulsion_strength': {'units': 'micron/min', 'value': 10.0}, 'relative_maximum_adhesion_distance': {'units': 'dimensionless', 'value': 1.25}}
			cell.dict['oxygen']={'secretion_rate': 0.0, 'secretion_unit': '1/min', 'secretion_target': 38.0, 'uptake_rate': 10.0, 'uptake_unit': '1/min', 'net_export': 0.0, 'net_export_unit': 'total substrate/min', 'secretion_rate_MCS': 0.0, 'secretion_comment': '#WARNING: PhysiCell has a concept of "target secretion" that CompuCell3D does not. \n#The translating program attempts to implement it, but it may not be a 1 to 1 conversion.', 'net_export_MCS': 0.0, 'net_secretion_comment': '', 'uptake_rate_MCS': 1.0, 'uptake_comment': '#WARNING: To avoid negative concentrations, in CompuCell3D uptake is "bounded." \n# If the amount that would be uptaken is larger than the value at that pixel,\n# the uptake will be a set ratio of the amount available.\n# The conversion program uses 1 as the ratio,\n# you may want to revisit this.'}
			cell.dict['immunostimulatory_factor']={'secretion_rate': 0.0, 'secretion_unit': '1/min', 'secretion_target': 1.0, 'uptake_rate': 0.0, 'uptake_unit': '1/min', 'net_export': 0.0, 'net_export_unit': 'total substrate/min', 'secretion_rate_MCS': 0.0, 'secretion_comment': '#WARNING: PhysiCell has a concept of "target secretion" that CompuCell3D does not. \n#The translating program attempts to implement it, but it may not be a 1 to 1 conversion.', 'net_export_MCS': 0.0, 'net_secretion_comment': '', 'uptake_rate_MCS': 0.0, 'uptake_comment': '#WARNING: To avoid negative concentrations, in CompuCell3D uptake is "bounded." \n# If the amount that would be uptaken is larger than the value at that pixel,\n# the uptake will be a set ratio of the amount available.\n# The conversion program uses 1 as the ratio,\n# you may want to revisit this.'}





class SecretionSteppable(SteppableBasePy):

	def __init__(self, frequency=1):
		SteppableBasePy.__init__(self,frequency)

	def start(self):
		"""
		Called before MCS=0 while building the initial simulation
		"""
		self.secretors = {'oxygen': self.get_field_secretor('oxygen'),'immunostimulatory_factor': self.get_field_secretor('immunostimulatory_factor')}


	def step(self, mcs):
		"""
		Called every frequency MCS while executing the simulation

		:param mcs: current Monte Carlo step
		"""


		for field_name, secretor in self.secretors.items():
			for cell in self.cell_list_by_type(self.CANCER_CELL):
				if field_name in cell.dict.keys():
					data=cell.dict[field_name]
					seen = secretor.amountSeenByCell(cell)
#WARNING: PhysiCell has a concept of "target secretion" that CompuCell3D does not. 
#The translating program attempts to implement it, but it may not be a 1 to 1 conversion.
					net_secretion = max(0, seen-data['secretion_rate_MCS']) + data['net_export_MCS']
					# In PhysiCell cells are point-like, in CC3D they have an arbitrary shape. With this 
					# CC3D allows several different secretion locations: over the whole cell (what the translator uses),
					# just inside the cell surface, just outside the surface, at the surface. You should explore the options
					secretor.secreteInsideCell(cell, net_secretion)
			for cell in self.cell_list_by_type(self.IMMUNE_CELL):
				if field_name in cell.dict.keys():
					data=cell.dict[field_name]
					seen = secretor.amountSeenByCell(cell)
#WARNING: PhysiCell has a concept of "target secretion" that CompuCell3D does not. 
#The translating program attempts to implement it, but it may not be a 1 to 1 conversion.
					net_secretion = max(0, seen-data['secretion_rate_MCS']) + data['net_export_MCS']
					# In PhysiCell cells are point-like, in CC3D they have an arbitrary shape. With this 
					# CC3D allows several different secretion locations: over the whole cell (what the translator uses),
					# just inside the cell surface, just outside the surface, at the surface. You should explore the options
					secretor.secreteInsideCell(cell, net_secretion)


	def finish(self):
		"""
		Called after the last MCS to wrap up the simulation. Good place to close files and do post-processing
		"""


	def on_stop(self):
		"""
		Called if the simulation is stopped before the last MCS
		"""
		self.finish()


'''
    names = get_steppables_names(s)
    print(names)
