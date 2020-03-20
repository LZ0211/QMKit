# coding=utf-8
from .Atom import Atom
from .Bond import Bond
import re

class Atoms(object):
    def __init__(self):
        self.atoms = []
        self.hash = {}

    def addAtom(self,atom):
        if atom == None:
            return
        self.hash[atom.label] = len(self.atoms)
        self.atoms.append(atom)
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
        return self

    def removeAtomByLabel(self,label):
        if label in self.hash:
            self.atoms.pop(self.hash[label])
        self.updateHash()
        return self

    def updateHash(self):
        self.hash.clear()
        for i in range(len(self.atoms)):
            self.hash[self.atoms[i].label] = i
        return self

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
