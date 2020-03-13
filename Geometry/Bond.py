# coding=utf-8
class Bond(object):
    table = {
    }
    def __init__(self,type=1,length=1):
        self.type = type
        self.length = length

    def getBondType(self):
        return self.type

    def setBondType(self,type):
        self.type = type
        return self

    def getBondLength(self):
        return self.length

    def setBondLength(self,length):
        self.length = length
        return self
