from .Atoms import Atoms

class Fragment(Atoms):
    def __init__(self):
        super().__init__()
        self.bonds = []
        self.activated = None

    def active(self,label):
        pass
