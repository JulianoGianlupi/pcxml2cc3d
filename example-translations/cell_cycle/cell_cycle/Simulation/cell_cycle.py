from cc3d import CompuCellSetup
#------

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
from cell_cycleSteppables import ConstraintsSteppable
CompuCellSetup.register_steppable(steppable=ConstraintsSteppable())

from cell_cycleSteppables import PhenotypeSteppable
CompuCellSetup.register_steppable(steppable=PhenotypeSteppable())


CompuCellSetup.run()
