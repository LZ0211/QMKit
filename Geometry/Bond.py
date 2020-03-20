# coding=utf-8
from .Atom import Atom

class Bond(object):
    table = {
    }
    def __init__(self,atomA,atomB,type=1):
        self.head = atomA
        self.tail = atomB
        self.length = Atom.length(self.head,self.tail)
        self.type = type

    def getBondType(self):
        return self.type

    def setBondType(self,type):
        self.type = type
        return self

    def setHeadAtom(self,atom):
        self.head = atom
        self.length = Atom.length(self.head,self.tail)

    def setTailAtom(self,atom):
        self.tail = atom
        self.length = Atom.length(self.head,self.tail)

