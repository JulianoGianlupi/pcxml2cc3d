
from cc3d import CompuCellSetup
        

from testSteppables import testSteppable

CompuCellSetup.register_steppable(steppable=testSteppable(frequency=1))


CompuCellSetup.run()
