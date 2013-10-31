#!/usr/bin/env python
# encoding: utf-8
"""
untitled.py

Created by Andrew Ning on 2013-10-30.
Copyright (c) NREL. All rights reserved.
"""


import numpy as np
from ctypes import POINTER, c_int, c_double, c_float, Structure, byref, cast

c_int_p = POINTER(c_int)
c_double_p = POINTER(c_double)


def ip(x):
    return x.ctypes.data_as(c_int_p)


def dp(x):
    return x.ctypes.data_as(c_double_p)



class Nodes(Structure):
    _fields_ = [('nN', c_int),
                ('N', c_int_p),
                ('x', c_double_p),
                ('y', c_double_p),
                ('z', c_double_p),
                ('r', c_double_p)]

    @classmethod
    def setup(cls, N, x, y, z, r):
        return cls(len(N), ip(N), dp(x), dp(y), dp(z), dp(r))



class Reactions(Structure):
    _fields_ = [('nR', c_int),
                ('N', c_int_p),
                ('Rx', c_int_p),
                ('Ry', c_int_p),
                ('Rz', c_int_p),
                ('Rxx', c_int_p),
                ('Ryy', c_int_p),
                ('Rzz', c_int_p)]

    @classmethod
    def setup(cls, N, Rx, Ry, Rz, Rxx, Ryy, Rzz):
        return cls(len(N), ip(N), ip(Rx), ip(Ry), ip(Rz),
                   ip(Rxx), ip(Ryy), ip(Rzz))


class Elements(Structure):
    _fields_ = [('nE', c_int),
                ('EL', c_int_p),
                ('N1', c_int_p),
                ('N2', c_int_p),
                ('Ax', c_double_p),
                ('Asy', c_double_p),
                ('Asz', c_double_p),
                ('Jx', c_double_p),
                ('Iy', c_double_p),
                ('Iz', c_double_p),
                ('E', c_double_p),
                ('G', c_double_p),
                ('roll', c_double_p),
                ('density', c_double_p)]

    @classmethod
    def setup(cls, EL, N1, N2, Ax, Asy, Asz, Jx, Iy, Iz, E, G, roll, density):
        return cls(len(EL), ip(EL), ip(N1), ip(N2), dp(Ax), dp(Asy), dp(Asz),
                   dp(Jx), dp(Iy), dp(Iz), dp(E), dp(G), dp(roll), dp(density))



class PointLoads(Structure):
    _fields_ = [('nF', c_int),
                ('N', c_int_p),
                ('Fx', c_double_p),
                ('Fy', c_double_p),
                ('Fz', c_double_p),
                ('Mxx', c_double_p),
                ('Myy', c_double_p),
                ('Mzz', c_double_p)]

    @classmethod
    def setup(cls, N, Fx, Fy, Fz, Mxx, Myy, Mzz):
        return cls(len(N), ip(N), dp(Fx), dp(Fy), dp(Fz),
                   dp(Mxx), dp(Myy), dp(Mzz))



class UniformLoads(Structure):
    _fields_ = [('nU', c_int),
                ('EL', c_int_p),
                ('Ux', c_double_p),
                ('Uy', c_double_p),
                ('Uz', c_double_p)]

    @classmethod
    def setup(cls, EL, Ux, Uy, Uz):
        return cls(len(EL), ip(EL), dp(Ux), dp(Uy), dp(Uz))


class PrescribedDisplacements(Structure):
    _fields_ = [('nD', c_int),
                ('N', c_int_p),
                ('Dx', c_double_p),
                ('Dy', c_double_p),
                ('Dz', c_double_p),
                ('Dxx', c_double_p),
                ('Dyy', c_double_p),
                ('Dzz', c_double_p)]

    @classmethod
    def setup(cls, N, Dx, Dy, Dz, Dxx, Dyy, Dzz):
        return cls(len(N), ip(N), dp(Dx), dp(Dy), dp(Dz),
                   dp(Dxx), dp(Dyy), dp(Dzz))



class LoadCase(Structure):
    _fields_ = [('gx', c_double),
                ('gy', c_double),
                ('gz', c_double),
                ('pointLoads', PointLoads),
                ('uniformLoads', UniformLoads),
                ('prescribedDisplacements', PrescribedDisplacements)]

    @classmethod
    def setup(cls, gx, gy, gz, pL, uL, pD):
        return cls(gx, gy, gz, pL, uL, pD)




class Displacements(Structure):
    _fields_ = [('node', c_int_p),
                ('x', c_double_p),
                ('y', c_double_p),
                ('z', c_double_p),
                ('xrot', c_double_p),
                ('yrot', c_double_p),
                ('zrot', c_double_p)]

    @classmethod
    def setup(cls, nodes, x, y, z, xrot, yrot, zrot):
        return cls(ip(nodes), dp(x), dp(y), dp(z),
                   dp(xrot), dp(yrot), dp(zrot))


class Forces(Structure):
    _fields_ = [('element', c_int_p),
                ('node', c_int_p),
                ('Nx', c_double_p),
                ('Vy', c_double_p),
                ('Vz', c_double_p),
                ('Txx', c_double_p),
                ('Myy', c_double_p),
                ('Mzz', c_double_p)]

    @classmethod
    def setup(cls, element, node, Nx, Vy, Vz, Txx, Myy, Mzz):
        return cls(ip(element), ip(node), dp(Nx), dp(Vy), dp(Vz),
                   dp(Txx), dp(Myy), dp(Mzz))



class ReactionForces(Structure):
    _fields_ = [('node', c_int_p),
                ('Fx', c_double_p),
                ('Fy', c_double_p),
                ('Fz', c_double_p),
                ('Mxx', c_double_p),
                ('Myy', c_double_p),
                ('Mzz', c_double_p)]

    @classmethod
    def setup(cls, node, Fx, Fy, Fz, Mxx, Myy, Mzz):
        return cls(ip(node), dp(Fx), dp(Fy), dp(Fz),
                   dp(Mxx), dp(Myy), dp(Mzz))



# nodes

N = np.arange(1, 13, dtype=np.int32)
x = np.array([0.0, 120.0, 240.0, 360.0, 480.0, 600.0, 720.0, 120.0, 240.0, 360.0, 480.0, 600.0])
y = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 120.0, 120.0, 120.0, 120.0, 120.0])
z = np.zeros(12)
r = np.zeros(12)

nodes = Nodes.setup(N, x, y, z, r)


# reactions
NR = np.arange(1, 13, dtype=np.int32)
Rx = np.array([1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0], dtype=np.int32)
Ry = np.array([1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0], dtype=np.int32)
Rz = np.ones(12, dtype=np.int32)
Rxx = np.ones(12, dtype=np.int32)
Ryy = np.ones(12, dtype=np.int32)
Rzz = np.zeros(12, dtype=np.int32)

reactions = Reactions.setup(NR, Rx, Ry, Rz, Rxx, Ryy, Rzz)


# elements
EL = np.arange(1, 22, dtype=np.int32)
N1 = np.array([1, 2, 3, 4, 5, 6, 1, 2, 2, 3, 4, 4, 4, 5, 6, 6, 7, 8, 9, 10, 11], dtype=np.int32)
N2 = np.array([2, 3, 4, 5, 6, 7, 8, 8, 9, 9, 9, 10, 11, 11, 11, 12, 12, 9, 10, 11, 12], dtype=np.int32)
Ax = 10.0*np.ones(21)
Asy = 1.0*np.ones(21)
Asz = 1.0*np.ones(21)
Jx = 1.0*np.ones(21)
Iy = 1.0*np.ones(21)
Iz = 0.01*np.ones(21)
E = 29000*np.ones(21)
G = 11500*np.ones(21)
roll = np.zeros(21)
density = 7.33e-7*np.ones(21)

elements = Elements.setup(EL, N1, N2, Ax, Asy, Asz, Jx, Iy, Iz, E, G, roll, density)


# parameters
shear = 0               # 1: include shear deformation
geom = 0                # 1: include geometric stiffness
exagg_static = 10.0     # exaggerate mesh deformations
dx = 10.0               # x-axis increment for internal forces


# load cases

gx = 0.0
gy = -386.4
gz = 0.0

nF = np.array([2, 3, 4, 5, 6], dtype=np.int32)
Fx = np.zeros(5)
Fy = np.array([-10.0, -20.0, -20.0, -10.0, -20.0])
Fz = np.zeros(5)
Mxx = np.zeros(5)
Myy = np.zeros(5)
Mzz = np.zeros(5)

pL = PointLoads.setup(nF, Fx, Fy, Fz, Mxx, Myy, Mzz)


ELU = np.array([])
Ux = np.array([])
Uy = np.array([])
Uz = np.array([])
uL = UniformLoads.setup(ELU, Ux, Uy, Uz)

nD = np.array([8], dtype=np.int32)
Dx = np.array([0.1])
Dy = np.array([0.0])
Dz = np.array([0.0])
Dxx = np.array([0.0])
Dyy = np.array([0.0])
Dzz = np.array([0.0])
pD = PrescribedDisplacements.setup(nD, Dx, Dy, Dz, Dxx, Dyy, Dzz)

l1 = LoadCase.setup(gx, gy, gz, pL, uL, pD)
nCases = 1

cases = (LoadCase * nCases)(l1)


# outputs
nN = len(N)
nE = len(EL)

ddnodes = np.zeros(nN, dtype=np.int32)
ddx = np.zeros(nN)
ddy = np.zeros(nN)
ddz = np.zeros(nN)
ddxrot = np.zeros(nN)
ddyrot = np.zeros(nN)
ddzrot = np.zeros(nN)

disp1 = Displacements.setup(ddnodes, ddx, ddy, ddz, ddxrot, ddyrot, ddzrot)



Eelement = np.zeros(nE*2, dtype=np.int32)
Enode = np.zeros(nE*2, dtype=np.int32)
ENx = np.zeros(nE*2)
EVy = np.zeros(nE*2)
EVz = np.zeros(nE*2)
ETxx = np.zeros(nE*2)
EMyy = np.zeros(nE*2)
EMzz = np.zeros(nE*2)

forces1 = Forces.setup(Eelement, Enode, ENx, EVy, EVz, ETxx, EMyy, EMzz)


Rnode = np.zeros(nN, dtype=np.int32)
RFx = np.zeros(nN)
RFy = np.zeros(nN)
RFz = np.zeros(nN)
RMxx = np.zeros(nN)
RMyy = np.zeros(nN)
RMzz = np.zeros(nN)

reactionForces1 = ReactionForces.setup(Rnode, RFx, RFy, RFz, RMxx, RMyy, RMzz)

disp = (Displacements * nCases)(disp1)
forces = (Forces * nCases)(forces1)
reactionForces = (ReactionForces * nCases)(reactionForces1)

_frame3dd = np.ctypeslib.load_library('_pyframe3dd', '.')

_frame3dd.readinputs.argtypes = [POINTER(Nodes), POINTER(Reactions), POINTER(Elements),
    c_int, c_int, c_double, c_float, c_int, POINTER(LoadCase),
    POINTER(Displacements), POINTER(Forces), POINTER(ReactionForces)]
_frame3dd.readinputs.restype = None


_frame3dd.readinputs(nodes, reactions, elements, shear, geom, exagg_static, dx, nCases, cases,
    disp, forces, reactionForces)

# x = np.ctypeslib.ndpointer(dtype=c_double, shape=(nN,))

# print disp[0].x
# x = np.ctypeslib.as_array(disp[0].x, shape=(nN,))

print ddx
print ddy
print ddz
print ddxrot
print ddyrot
print ddzrot

print

print Eelement
print Enode
print ENx
print EVy
print EVz
print ETxx
print EMyy
print EMzz

print

print Rnode
print RFx
print RFy
print RFz
print RMxx
print RMyy
print RMzz

# 0.0
# 0.010776
# 0.035528
# 0.060279
# 0.086295
# 0.112311
# 0.129754
# 0.100000
# 0.089226
# 0.059394
# 0.029563
# 0.012122



# import os
# import numpy as np
# from ctypes import *

# # open('tmp.c', 'w').write('''\
# # typedef struct {
# #     int v1; double *v2;
# # } darray;
# # int test(darray *input, darray *output) {
# #     int i;
# #     /* note: this should first test for compatible size */
# #     for (i=0; i < input->v1; i++)
# #         *(output->v2 + i) = *(input->v2 + i) * 2;
# #     return 0;
# # }
# # ''')
# # os.system('gcc -shared -o tmp.so tmp.c')


# c_double_p = POINTER(c_double)


# class LoadCase(Structure):
#     _fields_ = [('n', c_int),
#                 ('f', POINTER(c_double))]

#     @classmethod
#     def fromnp(cls, a):
#         return cls(len(a), a.ctypes.data_as(POINTER(c_double)))


# # lib = CDLL('./tmp.so')
# # lib = CDLL('./pyframe3dd.so')
# f = np.ctypeslib.load_library('pyframe3dd', '.')
# # lib.readinputs.argtypes = c_int, POINTER(LoadCase)
# f.readinputs.argtypes = [c_int, POINTER(LoadCase)]
# f.readinputs.restype = None

# a1 = np.arange(3) + 1.0
# a2 = np.zeros(3)
# print 'before:', '\na1 =', a1, '\na2 =', a2

# # lib.test(darray.fromnp(a1), darray.fromnp(a2))
# # f.readinputs(1, LoadCase.fromnp(a1))
# a1 = np.linspace(1.0, 2.0, 5)
# a2 = np.linspace(11.0, 12.0, 7)
# l1 = LoadCase.fromnp(a1)
# l2 = LoadCase.fromnp(a2)

# l = (LoadCase * 2)(l1, l2)

# f.readinputs(2, l)
# # print 'after:', '\na1 =', a1, '\na2 =', a2


# exit()



# from ctypes import CDLL, c_double, POINTER, c_int, Structure, pointer, byref, cast
# import numpy as np

# f = np.ctypeslib.load_library('pyframe3dd', '.')


# # dat = [[1126877361,'sunny'], [1126877371,'rain'], [1126877385,'damn nasty'], [1126877387,'sunny']]
# # dat_dtype = np.dtype([('timestamp','i4'),('desc','|S12')])
# # arr = np.rec.fromrecords(dat,dtype=dat_dtype)

# # print arr


# # exit()

# # dtype = np.dtype([('nF', np.int), ('fx', '5f8')])
# # print dtype


# def asArray(v):
#     return v.ctypes.data_as(POINTER(c_double))


# class LoadCase(Structure):
#     _fields_ = [('nF', c_int),
#                 ('fx', POINTER(c_double))]

#     @classmethod
#     def fromnp(cls, a):
#         return cls(len(a), a.ctypes.data_as(POINTER(c_double)))


# # # l1 = loads()
# # # l2 = loads()
# # # l[0] = l1

# # # l[0].nF = 5
# # # l[0].fx = asArray(np.linspace(1.0, 3.0, 5))
# # # l[1].nF = 7
# # # l[1].fx = asArray(np.linspace(2.0, 10.0, 7))
# # l0 = LoadCase(5, asArray(np.linspace(1.0, 3.0, 5)))
# # # l0.fx = cast(l0.fx, POINTER(c_double))
# # # l0.fx = asArray(np.linspace(1.0, 3.0, 5))
# # # l1 = loads(7, asArray(np.linspace(2.0, 10.0, 7)))
# # # l = (loads * 2)(l0, l1)
# # # l = cast(l, POINTER(loads))

# # # print np.linspace(1.0, 3.0, 5)

# # print l0.nF
# # for i in range(5):
# #     print l0.fx[i]

# # exit()
# # # l = cast(l, POINTER(loads))

# f.readinputs.argtypes = [c_int, POINTER(LoadCase)]
# # f.readinputs.argtypes = [c_int, np.ctypeslib.ndpointer(dtype)]
# f.readinputs.restype = None

# # x = np.array([(1.0, 2), (3.0, 4)], dtype=[('x', float), ('y', int)])


# # # # print l[0].fx, exit()
# # y = np.array([(1, 1.0), (1, 1.0), (1, 1.0)], dtype=[('foo', '<i4'), ('bar', '<f8')])

# # x = np.recarray([(1, 2.0), (2, 3.0)], dtype=[('n', np.int), ('x', np.float)])

# # l = np.array([(5, np.linspace(1.0, 3.0, 5))], dtype=dtype)

#     # , (7, np.linspace(10, 20.0, 7))], dtype=dtype)

# f.readinputs(1, LoadCase.fromnp(np.linspace(1.0, 2.0, 5)))


# # # f.readinputs.argtypes = [c_int, POINTER(c_double)]

# # # v = np.linspace(1.1, 2.3, 10)



# # # f.readinputs(len(v), asArray(v))
