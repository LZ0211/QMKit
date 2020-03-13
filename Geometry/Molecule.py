# coding=utf-8
from .Atom import Atom
import re

class Molecule(object):
    def __init__(self,name="Unnamed",charge=0,multiplicity=1):
        self.name = name
        self.charge = charge
        self.multiplicity = multiplicity
        self.atoms = []
        self.links = {}
        self.hash = {}
        self.cache = {}

    def addAtom(self,atom):
        if atom == None:
            return
        self.hash[atom.label] = len(self.atoms)
        self.atoms.append(atom)
        self.updateMultiplicity()
        self.clearCache()
        return self

    def addlink(self,a,b,value=1):
        if not a.label in self.links:
            self.links[a.label] = {}
        if not b.label in self.links:
            self.links[b.label] = {}
        self.links[a.label][b.label] = value
        self.links[b.label][a.label] = value
        return self

    def queryAtom(self,label):
        return self.atoms[self.hash[label]]

    def queryAtomsBySymbol(self,symbol):
        return list(filter(lambda x:x.symbol==symbol,self.atoms))

    def removeAtom(self,atom):
        for i in range(len(self.atoms)):
            if self.atoms[i] == atom:
                self.atoms.pop(i)
                break
        self.updateHash()
        self.updateLink()
        self.clearCache()
        return self

    def removeAtomByLabel(self,label):
        if label in self.hash:
            self.atoms.pop(self.hash[label])
        self.updateHash()
        self.updateLink()
        self.clearCache()
        return self

    def updateHash(self):
        self.hash.clear()
        for i in range(len(self.atoms)):
            self.hash[self.atoms[i].label] = i
        return self

    def updateLink(self):
        for (k,v) in self.links.items():
            if not k in self.hash:
                del self.links[k]
            else:
                for (_k,_v) in v.items():
                    if not _k in self.hash:
                        del v[_k]
        return self

    def updateMultiplicity(self):
        electrons = 0-self.charge
        for atom in self.atoms:
            electrons += atom.electrons
        self.multiplicity = 2 if electrons % 2 else 1
    
    def clearCache(self):
        self.cache = {}

    def setCharge(self,charge):
        self.charge = int(charge)
        self.updateMultiplicity()

    def setMultiplicity(self,multiplicity):
        multiplicity = int(multiplicity)
        if (multiplicity - self.multiplicity) % 2 == 0:
            self.multiplicity = multiplicity
        else:
            raise Exception('multiplicity of %s is impossible!' % multiplicity)

    def translate(self,vector=(0,0,0)):
        for atom in self.atoms:
            atom.translate(*vector)
        return self
    
    def rotate(self,axis,angle):
        for atom in self.atoms:
            atom.rotate(axis,angle,(0,0,0))
        return self

    @property
    def size(self):
        return len(self.atoms)

    @property
    def weight(self):
        if 'weight' in self.cache:
            return self.cache['weight']
        weight = 0
        for atom in self.atoms:
            weight += atom.weight
        self.cache['weight'] = weight
        return weight

    @property
    def electrons(self):
        if 'electrons' in self.cache:
            return self.cache['electrons']
        electrons = 0
        for atom in self.atoms:
            electrons += atom.electrons
        electrons -= self.charge
        self.cache['weight'] = electrons
        return electrons

    def string(self,format='xyz'):
        if format == 'xyz':
            return self.toXYZ()
        elif format == 'gjf':
            return self.toGJF()
        elif format == 'pdb':
            return self.toPDB()
        elif format == 'mol':
            return self.toMOL()
        elif format == 'mol2':
            return self.toMOL2()

    def toXYZ(self):
        xyz = '%s\n%s\n' % (len(self.atoms),self.name)
        for atom in self.atoms:
            xyz += atom.string() + '\n'
        return xyz

    def toGJF(self):
        gjf = '# hf/3-21g\n\n%s\n\n%s %s\n' % (self.name,self.charge,self.multiplicity)
        isDynamic = False
        for atom in self.atoms:
            if atom.frozen == True:
                isDynamic = True
                break
        if isDynamic == False:
            for atom in self.atoms:
                gjf += atom.string() + '\n'
        else:
            for atom in self.atoms:
                if atom.frozen:
                    gjf += atom.string(" {symbol:>2s}  -1  {x:>-22.15f} {y:>-22.15f} {z:>-22.15f}") + '\n'
                else:
                    gjf += atom.string(" {symbol:>2s}   0  {x:>-22.15f} {y:>-22.15f} {z:>-22.15f}") + '\n'
        return gjf

    def toPDB(self):
        pdb = 'TITLE      %s\n' % self.name
        for i in range(self.size):
            atom = self.atoms[i]
            atom.index = i+1
            pdb += atom.string('HETATM{index:>5} {symbol:>2}           0    {x:>-8.3f}{y:>-8.3f}{z:>-8.3f}                      {symbol:>2}\n')
        pdb += 'END\n'
        for (k,v) in self.links.items():
            pdb += 'CONECT %5d' % self.queryAtom(k).index
            for id in v.keys():
                pdb += '%5d' % self.queryAtom(id).index
            pdb += '\n'
        pdb += '\n'
        return pdb

    def toMOL(self):
        pass

    def toMOL2(self):
        pass

    @staticmethod
    def parse(str,format='xyz'):
        pass

    @staticmethod
    def parseXYZ(str):
        lines = str.split('\n')
        atomCount = int(lines[0])
        mol = Molecule(lines[1])
        for line in lines[2:2+atomCount]:
            symbol, x, y, z = line.split()[:4]
            atom = Atom(symbol,pos=[x,y,z])
            mol.addAtom(atom)
        return mol

    @staticmethod
    def parsePDB(str):
        pass
    
    @staticmethod
    def parseGRO(str):
        lines = str.split('\n')
        atomCount = int(lines[1])
        name = lines[2].split()[0][1:]
        mol = Molecule(name)
        for line in lines[2:2+atomCount]:
            label, idx, x, y, z = line.split()[1:]
            atom = Atom(label[0:1],pos=[10*float(x),10*float(y),10*float(z)])
            if atom == None:
                atom = Atom(label[0:2],pos=[10*float(x),10*float(y),10*float(z)])
            if atom == None:
                #raise Exception('Unkonow element of %s' % label)
                return
            atom.setLabel(label)
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
        mol = Molecule(lines[start])
        charge,multiplicity = lines[start+2].split()
        for line in lines[start+3:]:
            cols = line.split()[:5]
            if len(cols) == 5:
                symbol, fix, x, y, z = cols
                symbol = re.sub(r'\(.*?\)','',symbol)
                atom = Atom(symbol,pos=[x,y,z])
                if fix == '-1':
                    atom.freeze()
                mol.addAtom(atom)
            elif len(cols) == 4:
                symbol, x, y, z = cols
                atom = Atom(symbol,pos=[x,y,z])
                mol.addAtom(atom)
        mol.setCharge(charge)
        mol.setMultiplicity(multiplicity)
        return mol
