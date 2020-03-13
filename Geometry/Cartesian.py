# coding=utf-8
import numpy as np
import math

class Cartesian(object):
    def __init__(self,x=0,y=0,z=0):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

    def __call__(self):
        return self.pos()

    def __sub__(self,other):
        return self.pos() - other.pos()

    def pos(self):
        return np.array([self.x,self.y,self.z])
    
    def center(self):
        self.x = 0
        self.y = 0
        self.z = 0
        return self

    def transform(self,matrix):
        vector = self.pos()
        self.x = np.dot(vector,matrix[0][0:3]) + matrix[0][3]
        self.y = np.dot(vector,matrix[1][0:3]) + matrix[1][3]
        self.z = np.dot(vector,matrix[2][0:3]) + matrix[2][3]
        return self

    def translate(self,a,b,c):
        return self.transform([
            [1,0,0,a],
            [0,1,0,b],
            [0,0,1,c]
        ])

    def rotate(self,axis,angle,center=[0,0,0]):
        angle = angle / 180 * math.pi
        axis = np.array(axis)
        axis = axis / np.linalg.norm(axis)
        (x,y,z) = axis
        cos = math.cos(angle)
        _cos = 1-cos
        sin = math.sin(angle)
        matrix = [
            [x*x*_cos+cos,x*y*_cos+z*sin,x*z*_cos-y*sin,0],
            [x*y*_cos-z*sin,y*y*_cos+cos,y*z*_cos+x*sin,0],
            [x*z*_cos+y*sin,y*z*_cos-x*sin,z*z*_cos+cos,0],
        ]
        matrix[0][3] =center[0]*(1-matrix[0][0])-center[1]*matrix[0][1]-center[2]*matrix[0][2]
        matrix[1][3] =center[1]*(1-matrix[1][1])-center[0]*matrix[1][0]-center[2]*matrix[1][2]
        matrix[2][3] =center[2]*(1-matrix[2][2])-center[0]*matrix[2][0]-center[1]*matrix[2][1]
        self.transform(matrix)
        return self

    def align(self,axis,center=[0,0,0]):
        axis = np.array(axis)
        center = np.array(center)
        vector = self.pos()-center
        cos = np.dot(axis,vector) / np.linalg.norm(vector) / np.linalg.norm(axis)
        sin = math.sqrt(1-cos*cos)
        axis = np.cross(axis,vector)
        axis = axis / np.linalg.norm(axis)
        (x,y,z) = axis
        _cos = 1-cos
        matrix = [
            [x*x*_cos+cos,x*y*_cos+z*sin,x*z*_cos-y*sin,0],
            [x*y*_cos-z*sin,y*y*_cos+cos,y*z*_cos+x*sin,0],
            [x*z*_cos+y*sin,y*z*_cos-x*sin,z*z*_cos+cos,0]
        ]
        matrix[0][3] =center[0]*(1-matrix[0][0])-center[1]*matrix[0][1]-center[2]*matrix[0][2]
        matrix[1][3] =center[1]*(1-matrix[1][1])-center[0]*matrix[1][0]-center[2]*matrix[1][2]
        matrix[2][3] =center[2]*(1-matrix[2][2])-center[0]*matrix[2][0]-center[1]*matrix[2][1]
        self.transform(matrix)
        return self

    def rotateX(self,angle,center=[0,0,0]):
        angle = angle / 180 * math.pi
        cos = math.cos(angle)
        sin = math.sin(angle)
        return self.transform([
            [1,0,0,0],
            [0,cos,-sin,center[1]*(1-cos)+center[2]*sin],
            [0,sin,cos,center[2]*(1-cos)-center[1]*sin]
        ])

    def rotateY(self,angle,center=[0,0,0]):
        angle = angle / 180 * math.pi
        cos = math.cos(angle)
        sin = math.sin(angle)
        return self.transform([
            [cos,0,sin,center[0]*(1-cos)-center[2]*sin],
            [0,1,0,0],
            [-sin,0,cos,center[2]*(1-cos)+center[0]*sin]
        ])

    def rotateZ(self,angle,center=[0,0,0]):
        angle = angle / 180 * math.pi
        cos = math.cos(angle)
        sin = math.sin(angle)
        return self.transform([
            [cos,-sin,0,center[0]*(1-cos)+center[1]*sin],
            [sin,cos,0,center[1]*(1-cos)-center[0]*sin],
            [0,0,1,0]
        ])

    def alignZ(self,center=[0,0,0]):
        center = np.array(center)
        vector = self.pos()-center
        (a,b,c) = vector
        r1 = math.sqrt(b*b+c*c)
        r2 = math.sqrt(a*a+b*b+c*c)
        matrix = [
            [r1/r2, -a*b/r1/r2, -a*c/r1/r2, 0],
            [0,     c/r1,       -b/r1,      0],
            [a/r2,  b/r2,       c/r2,       0]
        ]
        matrix[0][3] =center[0]*(1-matrix[0][0])-center[1]*matrix[0][1]-center[2]*matrix[0][2]
        matrix[1][3] =center[1]*(1-matrix[1][1])-center[0]*matrix[1][0]-center[2]*matrix[1][2]
        matrix[2][3] =center[2]*(1-matrix[2][2])-center[0]*matrix[2][0]-center[1]*matrix[2][1]
        return self.transform(matrix)

    def alignY(self,center=[0,0,0]):
        center = np.array(center)
        vector = self.pos()-center
        (a,b,c) = vector
        r1 = math.sqrt(b*b+c*c)
        r2 = math.sqrt(a*a+b*b+c*c)
        matrix = [
            [r1/r2, -a*b/r1/r2, -a*c/r1/r2, 0],
            [a/r2,  b/r2,       c/r2,       0],
            [0,     c/r1,       -b/r1,      0]
        ]
        matrix[0][3] =center[0]*(1-matrix[0][0])-center[1]*matrix[0][1]-center[2]*matrix[0][2]
        matrix[1][3] =center[1]*(1-matrix[1][1])-center[0]*matrix[1][0]-center[2]*matrix[1][2]
        matrix[2][3] =center[2]*(1-matrix[2][2])-center[0]*matrix[2][0]-center[1]*matrix[2][1]
        return self.transform(matrix)

    def alignX(self,center=[0,0,0]):
        center = np.array(center)
        vector = self.pos()-center
        (a,b,c) = vector
        r1 = math.sqrt(b*b+c*c)
        r2 = math.sqrt(a*a+b*b+c*c)
        matrix = [
            [a/r2,  b/r2,       c/r2,       0],
            [0,     c/r1,       -b/r1,      0],
            [r1/r2, -a*b/r1/r2, -a*c/r1/r2, 0]
        ]
        matrix[0][3] =center[0]*(1-matrix[0][0])-center[1]*matrix[0][1]-center[2]*matrix[0][2]
        matrix[1][3] =center[1]*(1-matrix[1][1])-center[0]*matrix[1][0]-center[2]*matrix[1][2]
        matrix[2][3] =center[2]*(1-matrix[2][2])-center[0]*matrix[2][0]-center[1]*matrix[2][1]
        return self.transform(matrix)

    @staticmethod
    def normal(v):
        return np.linalg.norm(v)

    @staticmethod
    def length(a,b):
        return np.linalg.norm(a-b)

    @staticmethod
    def angle(a,b,c):
        v1 = a-b
        v2 = c-b
        return math.acos(np.dot(v1,v2)/np.linalg.norm(v1)/np.linalg.norm(v2))/math.pi*180

    @staticmethod
    def dihedral(a,b,c,d):
        v1 = b-a
        v2 = c-b
        v3 = d-c
        a1 = np.cross(v1,v2)
        a2 = np.cross(v2,v3)
        b = -1
        if np.dot(a1,v3) >= 0:
            b = 1
        return b* math.acos(np.dot(a1,a2)/np.linalg.norm(a1)/np.linalg.norm(a2))/math.pi*180
