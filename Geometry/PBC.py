from .Molecule import Molecule
from .Cartesian import Cartesian
from .Atom import Atom
import numpy as np

class PBC(Molecule):
    def __init__(self,name="Unnamed",charge=0,multiplicity=1,vectors=[(10,0,0),(0,10,0),(0,0,10)]):
        super().__init__(name,charge,multiplicity)
        self.dimensions = len(vectors)
        self.vectors = vectors

    def setVectors(self,vectors):
        self.dimensions = len(vectors)
        self.vectors = vectors
        return self
        
    def setVectorX(self,vector):
        self.vectors[0] = vector
        return self

    def setVectorY(self,vector):
        if self.dimensions >= 2:
            self.vectors[1] = vector
        return self

    def setVectorZ(self,vector):
        if self.dimensions >=3:
            self.vectors[2] = vector
        return self

    def toGJF(self):
        gjf = Molecule.toGJF(self)
        for vector in self.vectors:
            vector = Cartesian(*vector)
            vector.symbol = 'Tv'
            gjf += Atom.string(vector) + '\n'
        return gjf

    def toVASP(self):
        hash = {}
        for atom in self.atoms:
            symbol = atom.symbol
            if not symbol in hash:
                hash[symbol] = []
            hash[symbol].append(atom)
        isDynamic = False
        for atom in self.atoms:
            if atom.frozen == True:
                isDynamic = True
                break
        vasp = '%s\n%22.15f\n' % (self.name,1)
        for vector in self.vectors:
            vasp += '%22.15f %22.15f %22.15f\n' % tuple(vector)
        for k in hash.keys():
            vasp += '%8s' % k
        vasp += '\n'
        for v in hash.values():
            vasp += '%8d' % len(v)
        vasp += '\n'
        if isDynamic:
            vasp += 'Selective dynamics\n'
        vasp += 'Cartesian\n'
        if isDynamic:
            for atoms in hash.values():
                for atom in atoms:
                    if atom.frozen:
                        vasp += atom.string("{x:>-22.15f} {y:>-22.15f} {z:>-22.15f}   F   F   F\n")
                    else:
                        vasp += atom.string("{x:>-22.15f} {y:>-22.15f} {z:>-22.15f}   T   T   T\n")
        else:
            for atoms in hash.values():
                for atom in atoms:
                    vasp += atom.string("{x:>-22.15f} {y:>-22.15f} {z:>-22.15f}\n")
        return vasp
    
    @staticmethod
    def parseXYZ(str):
        lines = str.split('\n')
        atomCount = int(lines[0])
        mol = PBC(lines[1])
        for line in lines[2:2+atomCount]:
            symbol, x, y, z = line.split()[:4]
            atom = Atom(symbol,pos=[x,y,z])
            mol.addAtom(atom)
        return mol

    @staticmethod
    def parseVASP(str):
        lines = str.split('\n')
        name = lines[0]
        mol = PBC(name)
        scale = float(lines[1])
        vectors = [lines[2].split(),lines[3].split(),lines[4].split()]
        for vector in vectors:
            vector[0] = scale*float(vector[0])
            vector[1] = scale*float(vector[1])
            vector[2] = scale*float(vector[2])
        mol.setVectors(vectors)
        matrix = np.array(vectors)
        elements = lines[5].split()
        counts = lines[6].split()
        atoms = []
        for i in range(len(elements)):
            count = counts[i]
            element = elements[i]
            atoms.extend([element for i in range(int(count))])
        start = 7
        isDynamic = ~lines[7].find('Selective')
        if isDynamic:
            start = 8
        isCartesian = (lines[start][0]=='C')
        start += 1
        for i in range(len(atoms)):
            line = lines[start+i]
            arr = line.split()
            pos = list(map(float,arr[:3]))
            if not isCartesian:
                pos = np.dot(matrix,pos)
            atom = Atom(atoms[i],pos=pos)
            if isDynamic and ~line.find('F'):
                atom.freeze()
            mol.addAtom(atom)
        return mol

    @staticmethod
    def parseGJF(str):
        lines = str.split('\n')
        start = 0
        for i in range(10):
            if lines[i] == '' and lines[i+2] == '':
                start = i+1
                break
        mol = PBC(lines[start])
        charge,multiplicity = lines[start+2].split()
        vectors = []
        for line in lines[start+3:]:
            cols = line.split()[:5]
            if len(cols) == 0:
                break
            if len(cols) == 5:
                symbol, fix, x, y, z = cols
                atom = Atom(symbol,pos=[x,y,z])
                if fix == '-1':
                    atom.freeze()
                mol.addAtom(atom)
            elif len(cols) == 4:
                symbol, x, y, z = cols
                if symbol.capitalize() == 'Tv':
                    vectors.append((float(x),float(y),float(z)))
                else:
                    atom = Atom(symbol,pos=[x,y,z])
                    mol.addAtom(atom)
        mol.setVectors(vectors)
        mol.setCharge(charge)
        mol.setMultiplicity(multiplicity)
        return mol

