from pathlib import Path
from os.path import join, isdir


def generate_steppable_file(path, steppable_fname, steppable_string):
    path = Path(path)

    if not isdir(path):
        path.mkdir(parents=True)

    with open(join(path, steppable_fname), "w+") as f:
        f.write(steppable_string.replace("\t", "    "))


if __name__ == "__main__":
    steppable_string = """
from cc3d.cpp.PlayerPython import *
from cc3d import CompuCellSetup
from cc3d.core.PySteppables import *
import numpy as np

class ConstraintsSteppable(SteppableBasePy):

	def __init__(self, frequency=1):
		SteppableBasePy.__init__(self,frequency)

	def start(self):
		return
		"""

    generate_steppable_file("D:/test_pc2cc3d/Simulation", "steppables.py", steppable_string)
