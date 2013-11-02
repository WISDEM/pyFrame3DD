#!/usr/bin/env python
# encoding: utf-8
"""
frame.py

Created by Andrew Ning on 2013-11-01.
Copyright (c) NREL. All rights reserved.
"""

import numpy as np
import math
from ctypes import POINTER, c_int, c_double, Structure


c_int_p = POINTER(c_int)
c_double_p = POINTER(c_double)


def ip(x):
    return x.ctypes.data_as(c_int_p)


def dp(x):
    return x.ctypes.data_as(c_double_p)



# --------------
# General Inputs
# --------------

class Nodes(Structure):
    _fields_ = [('nN', c_int),
                ('N', c_int_p),
                ('x', c_double_p),
                ('y', c_double_p),
                ('z', c_double_p),
                ('r', c_double_p)]



class Reactions(Structure):
    _fields_ = [('nR', c_int),
                ('N', c_int_p),
                ('Rx', c_int_p),
                ('Ry', c_int_p),
                ('Rz', c_int_p),
                ('Rxx', c_int_p),
                ('Ryy', c_int_p),
                ('Rzz', c_int_p)]



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



class OtherElementData(Structure):
    _fields_ = [('shear', c_int),
                ('geom', c_int),
                ('exagg_static', c_double),
                ('dx', c_double)]



# --------------
# Load Inputs
# --------------


class PointLoads(Structure):
    _fields_ = [('nF', c_int),
                ('N', c_int_p),
                ('Fx', c_double_p),
                ('Fy', c_double_p),
                ('Fz', c_double_p),
                ('Mxx', c_double_p),
                ('Myy', c_double_p),
                ('Mzz', c_double_p)]




class UniformLoads(Structure):
    _fields_ = [('nU', c_int),
                ('EL', c_int_p),
                ('Ux', c_double_p),
                ('Uy', c_double_p),
                ('Uz', c_double_p)]



class TrapezoidalLoads(Structure):
    _fields_ = [('nW', c_int),
                ('EL', c_int_p),
                ('xx1', c_double_p),
                ('xx2', c_double_p),
                ('wx1', c_double_p),
                ('wx2', c_double_p),
                ('xy1', c_double_p),
                ('xy2', c_double_p),
                ('wy1', c_double_p),
                ('wy2', c_double_p),
                ('xz1', c_double_p),
                ('xz2', c_double_p),
                ('wz1', c_double_p),
                ('wz2', c_double_p)]



class ElementLoads(Structure):
    _fields_ = [('nP', c_int),
                ('EL', c_int_p),
                ('Px', c_double_p),
                ('Py', c_double_p),
                ('Pz', c_double_p),
                ('x', c_double_p)]


class TemperatureLoads(Structure):
    _fields_ = [('nT', c_int),
                ('EL', c_int_p),
                ('a', c_double_p),
                ('hy', c_double_p),
                ('hz', c_double_p),
                ('Typ', c_double_p),
                ('Tym', c_double_p),
                ('Tzp', c_double_p),
                ('Tzm', c_double_p),
                ]




class PrescribedDisplacements(Structure):
    _fields_ = [('nD', c_int),
                ('N', c_int_p),
                ('Dx', c_double_p),
                ('Dy', c_double_p),
                ('Dz', c_double_p),
                ('Dxx', c_double_p),
                ('Dyy', c_double_p),
                ('Dzz', c_double_p)]




class LoadCase(Structure):
    _fields_ = [('gx', c_double),
                ('gy', c_double),
                ('gz', c_double),
                ('pointLoads', PointLoads),
                ('uniformLoads', UniformLoads),
                ('trapezoidalLoads', TrapezoidalLoads),
                ('elementLoads', ElementLoads),
                ('temperatureLoads', TemperatureLoads),
                ('prescribedDisplacements', PrescribedDisplacements)]


# --------------
# Dynamic Inputs
# --------------


class DynamicData(Structure):
    _fields_ = [('nM', c_int),
                ('Mmethod', c_int),
                ('lump', c_int),
                ('tol', c_double),
                ('shift', c_double),
                ('exagg_modal', c_double)]


class ExtraInertia(Structure):
    _fields_ = [('nI', c_int),
                ('N', c_int_p),
                ('EMs', c_double_p),
                ('EMx', c_double_p),
                ('EMy', c_double_p),
                ('EMz', c_double_p)]


class ExtraMass(Structure):
    _fields_ = [('nX', c_int),
                ('EL', c_int_p),
                ('EMs', c_double_p)]


class Condensation(Structure):
    _fields_ = [('Cmethod', c_int),
                ('nC', c_int),
                ('N', c_int_p),
                ('cx', c_double_p),
                ('cy', c_double_p),
                ('cz', c_double_p),
                ('cxx', c_double_p),
                ('cyy', c_double_p),
                ('czz', c_double_p),
                ('m', c_int_p)]



# --------------
# Static Data Outputs
# --------------

class Displacements(Structure):
    _fields_ = [('node', c_int_p),
                ('x', c_double_p),
                ('y', c_double_p),
                ('z', c_double_p),
                ('xrot', c_double_p),
                ('yrot', c_double_p),
                ('zrot', c_double_p)]



class Forces(Structure):
    _fields_ = [('element', c_int_p),
                ('node', c_int_p),
                ('Nx', c_double_p),
                ('Vy', c_double_p),
                ('Vz', c_double_p),
                ('Txx', c_double_p),
                ('Myy', c_double_p),
                ('Mzz', c_double_p)]




class ReactionForces(Structure):
    _fields_ = [('node', c_int_p),
                ('Fx', c_double_p),
                ('Fy', c_double_p),
                ('Fz', c_double_p),
                ('Mxx', c_double_p),
                ('Myy', c_double_p),
                ('Mzz', c_double_p)]




# --------------
# Internal Force Outputs
# --------------


class InternalForces(Structure):
    _fields_ = [('x', c_double_p),
                ('Nx', c_double_p),
                ('Vy', c_double_p),
                ('Vz', c_double_p),
                ('Tx', c_double_p),
                ('My', c_double_p),
                ('Mz', c_double_p),
                ('Dx', c_double_p),
                ('Dy', c_double_p),
                ('Dz', c_double_p),
                ('Rx', c_double_p),
                ]





# --------------
# Modal Outputs
# --------------


class MassResults(Structure):
    _fields_ = [('total_mass', c_double),
                ('struct_mass', c_double),
                ('N', c_int_p),
                ('xmass', c_double_p),
                ('ymass', c_double_p),
                ('zmass', c_double_p),
                ('xinrta', c_double_p),
                ('yinrta', c_double_p),
                ('zinrta', c_double_p),
                ]



class ModalResults(Structure):
    _fields_ = [('freq', c_double),
                ('xmpf', c_double),
                ('ympf', c_double),
                ('zmpf', c_double),
                ('N', c_int_p),
                ('xdsp', c_double_p),
                ('ydsp', c_double_p),
                ('zdsp', c_double_p),
                ('xrot', c_double_p),
                ('yrot', c_double_p),
                ('zrot', c_double_p),
                ]








class Frame(object):
    """docstring for Frame"""


    def __init__(self, N, x, y, z, r,
            NR, Rx, Ry, Rz, Rxx, Ryy, Rzz,
            EL, N1, N2, Ax, Asy, Asz, Jx, Iy, Iz, E, G, roll, density,
            shear, geom, exagg_static, dx):


        # copy to prevent user error with name clashes (c is sharing the address space)

        self.N = N.astype(np.int32)
        self.x = np.copy(x)
        self.y = np.copy(y)
        self.z = np.copy(z)
        self.r = np.copy(r)
        self.NR = NR.astype(np.int32)
        self.Rx = Rx.astype(np.int32)
        self.Ry = Ry.astype(np.int32)
        self.Rz = Rz.astype(np.int32)
        self.Rxx = Rxx.astype(np.int32)
        self.Ryy = Ryy.astype(np.int32)
        self.Rzz = Rzz.astype(np.int32)
        self.EL = EL.astype(np.int32)
        self.N1 = N1.astype(np.int32)
        self.N2 = N2.astype(np.int32)
        self.Ax = np.copy(Ax)
        self.Asy = np.copy(Asy)
        self.Asz = np.copy(Asz)
        self.Jx = np.copy(Jx)
        self.Iy = np.copy(Iy)
        self.Iz = np.copy(Iz)
        self.E = np.copy(E)
        self.G = np.copy(G)
        self.roll = np.copy(roll)
        self.density = np.copy(density)

        self.dx = dx


        # setup general data
        self.nodes = Nodes(len(N), ip(self.N), dp(self.x), dp(self.y),
            dp(self.z), dp(self.r))

        self.reactions = Reactions(len(NR), ip(self.NR), ip(self.Rx), ip(self.Ry),
            ip(self.Rz), ip(self.Rxx), ip(self.Ryy), ip(self.Rzz))

        self.elements = Elements(len(EL), ip(self.EL), ip(self.N1), ip(self.N2),
            dp(self.Ax), dp(self.Asy), dp(self.Asz), dp(self.Jx), dp(self.Iy),
            dp(self.Iz), dp(self.E), dp(self.G), dp(self.roll), dp(self.density))

        self.other = OtherElementData(shear, geom, exagg_static, dx)


        # create list for load cases
        self.loadCases = []
        self.staticOutputs = []


        # initialize no dynamics
        dynamic = DynamicAnalysis(nM=0, Mmethod=1, lump=0, tol=0.0, shift=0.0, exagg_modal=0.0)
        self.useDynamicAnalysis(dynamic)


        # load c module
        self._frame3dd = np.ctypeslib.load_library('_pyframe3dd', '.')

        self._frame3dd.run.argtypes = [POINTER(Nodes), POINTER(Reactions), POINTER(Elements),
            POINTER(OtherElementData), c_int, POINTER(LoadCase),
            POINTER(DynamicData), POINTER(ExtraInertia), POINTER(ExtraMass),
            POINTER(Condensation),
            POINTER(Displacements), POINTER(Forces), POINTER(ReactionForces),
            POINTER(POINTER(InternalForces)), POINTER(MassResults), POINTER(ModalResults)]

        self._frame3dd.run.restype = c_int





    def addLoadCase(self, loadCase):

        self.loadCases.append(loadCase)

        so = StaticOutputs(len(self.N), len(self.EL), self.x, self.y, self.z, self.N1, self.N2, self.dx)
        self.staticOutputs.append(so)


    def useDynamicAnalysis(self, dynamic):

        self.dynamic = dynamic
        self.dynamicOutputs = DynamicOutputs(len(self.N), dynamic.nM)


    def run(self):

        # load cases
        nCases = len(self.loadCases)

        noloads = False
        if nCases == 0:
            self.loadCases.append(StaticLoadCase(0.0, 0.0, 0.0))
            noloads = True

        loadcases = (LoadCase * nCases)()
        disp = (Displacements * nCases)()
        forces = (Forces * nCases)()
        reactions = (ReactionForces * nCases)()
        internalForces = (POINTER(InternalForces) * nCases)()

        for i in range(nCases):
            lci = self.loadCases[i]
            loadcases[i] = LoadCase(lci.gx, lci.gy, lci.gz, lci.pL,
                lci.uL, lci.tL, lci.eL, lci.tempL, lci.pD)
            disp[i] = self.staticOutputs[i].disp
            forces[i] = self.staticOutputs[i].forces
            reactions[i] = self.staticOutputs[i].reactions
            internalForces[i] = self.staticOutputs[i].internalForces

        d = self.dynamic
        do = self.dynamicOutputs



        self._frame3dd.run(self.nodes, self.reactions, self.elements, self.other,
            nCases, loadcases,
            d.dynamicData, d.extraInertia, d.extraMass, d.condensation,
            disp, forces, reactions, internalForces, do.massResults, do.modalResults)


        if noloads:
            self.loadCases = []


        so = self.staticOutputs[0]
        print so.dx
        print so.dy
        print so.dz
        print so.dxrot
        print so.dyrot
        print so.dzrot


class StaticLoadCase(object):
    """docstring"""


    def __init__(self, gx, gy, gz):

        self.gx = gx
        self.gy = gy
        self.gz = gz

        i = np.array([], dtype=np.int32)
        d = np.array([])

        self.changePointLoads(i, d, d, d, d, d, d)

        self.changeUniformLoads(i, d, d, d)

        self.changeTrapezoidalLoads(i, d, d, d, d, d, d, d, d, d, d, d, d)

        self.changeElementLoads(i, d, d, d, d)

        self.changeTemperatureLoads(i, d, d, d, d, d, d, d)

        self.changePrescribedDisplacements(i, d, d, d, d, d, d)




    def changePointLoads(self, N, Fx, Fy, Fz, Mxx, Myy, Mzz):

        self.NF = N.astype(np.int32)
        self.Fx = np.copy(Fx)
        self.Fy = np.copy(Fy)
        self.Fz = np.copy(Fz)
        self.Mxx = np.copy(Mxx)
        self.Myy = np.copy(Myy)
        self.Mzz = np.copy(Mzz)

        self.pL = PointLoads(len(N), ip(self.NF), dp(self.Fx), dp(self.Fy), dp(self.Fz),
                   dp(self.Mxx), dp(self.Myy), dp(self.Mzz))


    def changeUniformLoads(self, EL, Ux, Uy, Uz):

        self.ELU = EL.astype(np.int32)
        self.Ux = np.copy(Ux)
        self.Uy = np.copy(Uy)
        self.Uz = np.copy(Uz)

        self.uL = UniformLoads(len(EL), ip(self.ELU), dp(self.Ux), dp(self.Uy), dp(self.Uz))



    def changeTrapezoidalLoads(self, EL, xx1, xx2, wx1, wx2, xy1, xy2, wy1, wy2, xz1, xz2, wz1, wz2):

        self.ELT = EL.astype(np.int32)
        self.xx1 = np.copy(xx1)
        self.xx2 = np.copy(xx2)
        self.wx1 = np.copy(wx1)
        self.wx2 = np.copy(wx2)
        self.xy1 = np.copy(xy1)
        self.xy2 = np.copy(xy2)
        self.wy1 = np.copy(wy1)
        self.wy2 = np.copy(wy2)
        self.xz1 = np.copy(xz1)
        self.xz2 = np.copy(xz2)
        self.wz1 = np.copy(wz1)
        self.wz2 = np.copy(wz2)

        self.tL = TrapezoidalLoads(len(EL), ip(self.ELT), dp(self.xx1), dp(self.xx2), dp(self.wx1), dp(self.wx2),
            dp(self.xy1), dp(self.xy2), dp(self.wy1), dp(self.wy2), dp(self.xz1), dp(self.xz2), dp(self.wz1), dp(self.wz2))



    def changeElementLoads(self, EL, Px, Py, Pz, x):

        self.ELE = EL.astype(np.int32)
        self.Px = np.copy(Px)
        self.Py = np.copy(Py)
        self.Pz = np.copy(Pz)
        self.xE = np.copy(x)


        self.eL = ElementLoads(len(EL), ip(self.ELE), dp(self.Px), dp(self.Py), dp(self.Pz), dp(self.xE))



    def changeTemperatureLoads(self, EL, a, hy, hz, Typ, Tym, Tzp, Tzm):

        self.ELTemp = EL.astype(np.int32)
        self.a = np.copy(a)
        self.hy = np.copy(hy)
        self.hz = np.copy(hz)
        self.Typ = np.copy(Typ)
        self.Tym = np.copy(Tym)
        self.Tzp = np.copy(Tzp)
        self.Tzm = np.copy(Tzm)


        self.tempL = TemperatureLoads(len(EL), ip(self.ELTemp), dp(self.a), dp(self.hy), dp(self.hz), dp(self.Typ),
            dp(self.Tym), dp(self.Tzp), dp(self.Tzm))



    def changePrescribedDisplacements(self, N, Dx, Dy, Dz, Dxx, Dyy, Dzz):

        self.ND = N.astype(np.int32)
        self.Dx = np.copy(Dx)
        self.Dy = np.copy(Dy)
        self.Dz = np.copy(Dz)
        self.Dxx = np.copy(Dxx)
        self.Dyy = np.copy(Dyy)
        self.Dzz = np.copy(Dzz)

        self.pD = PrescribedDisplacements(len(N), ip(self.ND), dp(self.Dx), dp(self.Dy), dp(self.Dz),
                   dp(self.Dxx), dp(self.Dyy), dp(self.Dzz))




class DynamicAnalysis(object):
    """docstring"""


    def __init__(self, nM, Mmethod, lump, tol, shift, exagg_modal):

        self.nM = nM

        self.dynamicData = DynamicData(nM, Mmethod, lump, tol, shift, exagg_modal)

        i = np.array([], dtype=np.int32)
        d = np.array([])

        self.changeExtraInertia(i, d, d, d, d)

        self.changeExtraMass(i, d)

        self.changeCondensationData(0, i, d, d, d, d, d, d, i)



    def changeExtraInertia(self, N, EMs, EMx, EMy, EMz):

        self.NI = N.astype(np.int32)
        self.EMs = np.copy(EMs)
        self.EMx = np.copy(EMx)
        self.EMy = np.copy(EMy)
        self.EMz = np.copy(EMz)


        self.extraInertia = ExtraInertia(len(N), ip(self.NI), dp(self.EMs),
            dp(self.EMx), dp(self.EMy), dp(self.EMz))



    def changeExtraMass(self, EL, EMs):

        self.ELM = EL.astype(np.int32)
        self.EMsM = np.copy(EMs)

        self.extraMass = ExtraMass(len(EL), ip(self.ELM), dp(self.EMs))




    def changeCondensationData(self, Cmethod, N, cx, cy, cz, cxx, cyy, czz, m):

        self.NC = N.astype(np.int32)
        self.cx = np.copy(cx)
        self.cy = np.copy(cy)
        self.cz = np.copy(cx)
        self.cxx = np.copy(cxx)
        self.cyy = np.copy(cyy)
        self.czz = np.copy(czz)
        self.mC = m.astype(np.int32)

        self.condensation = Condensation(Cmethod, len(N), ip(self.NC), dp(self.cx), dp(self.cy), dp(self.cz),
            dp(self.cxx), dp(self.cyy), dp(self.czz), ip(self.mC))








class StaticOutputs(object):
    """docstring"""


    def __init__(self, nN, nE, x, y, z, N1, N2, dx):


        self.ND = np.zeros(nN, dtype=np.int32)
        self.dx = np.zeros(nN)
        self.dy = np.zeros(nN)
        self.dz = np.zeros(nN)
        self.dxrot = np.zeros(nN)
        self.dyrot = np.zeros(nN)
        self.dzrot = np.zeros(nN)

        self.disp = Displacements(ip(self.ND), dp(self.dx), dp(self.dy), dp(self.dz),
                   dp(self.dxrot), dp(self.dyrot), dp(self.dzrot))


        self.ELF = np.zeros(nE*2, dtype=np.int32)
        self.NF = np.zeros(nE*2, dtype=np.int32)
        self.Nx = np.zeros(nE*2)
        self.Vy = np.zeros(nE*2)
        self.Vz = np.zeros(nE*2)
        self.Txx = np.zeros(nE*2)
        self.Myy = np.zeros(nE*2)
        self.Mzz = np.zeros(nE*2)

        self.forces = Forces(ip(self.ELF), ip(self.NF), dp(self.Nx), dp(self.Vy), dp(self.Vz),
                   dp(self.Txx), dp(self.Myy), dp(self.Mzz))


        self.NR = np.zeros(nN, dtype=np.int32)
        self.RFx = np.zeros(nN)
        self.RFy = np.zeros(nN)
        self.RFz = np.zeros(nN)
        self.RMxx = np.zeros(nN)
        self.RMyy = np.zeros(nN)
        self.RMzz = np.zeros(nN)

        self.reactions = ReactionForces(ip(self.NR), dp(self.RFx), dp(self.RFy), dp(self.RFz),
                   dp(self.RMxx), dp(self.RMyy), dp(self.RMzz))



        self.fobj = [0]*nE
        self.internalForces = (InternalForces * nE)()

        for i in range(nE):

            self.fobj[i] = IForce(i, x, y, z, N1, N2, dx)
            self.internalForces[i] = self.fobj[i].iforce





class IForce(object):


    def __init__(self, i, x, y, z, N1, N2, dx):

        L = math.sqrt(
            (x[N2[i]-1] - x[N1[i]-1])**2 +
            (y[N2[i]-1] - y[N1[i]-1])**2 +
            (z[N2[i]-1] - z[N1[i]-1])**2
            )

        length = int(max(math.floor(L/dx), 1))

        self.x = np.zeros(length)
        self.Nx = np.zeros(length)
        self.Vy = np.zeros(length)
        self.Vz = np.zeros(length)
        self.Tx = np.zeros(length)
        self.My = np.zeros(length)
        self.Mz = np.zeros(length)
        self.Dx = np.zeros(length)
        self.Dy = np.zeros(length)
        self.Dz = np.zeros(length)
        self.Rx = np.zeros(length)

        self.iforce = InternalForces(dp(self.x), dp(self.Nx), dp(self.Vy), dp(self.Vz),
            dp(self.Tx), dp(self.My), dp(self.Mz), dp(self.Dx), dp(self.Dy), dp(self.Dz), dp(self.Rx))





class DynamicOutputs(object):
    """docstring"""


    def __init__(self, nN, nM):


        self.total_mass = 0.0
        self.struct_mass = 0.0

        self.NM = np.zeros(nN, dtype=np.int32)
        self.xmass = np.zeros(nN)
        self.ymass = np.zeros(nN)
        self.zmass = np.zeros(nN)
        self.xinrta = np.zeros(nN)
        self.yinrta = np.zeros(nN)
        self.zinrta = np.zeros(nN)

        self.massResults = MassResults(self.total_mass, self.struct_mass, ip(self.NM),
            dp(self.xmass), dp(self.ymass), dp(self.zmass),
            dp(self.xinrta), dp(self.yinrta), dp(self.zinrta))


        self.fobj = [0]*nM
        self.modalResults = (ModalResults * nM)()

        for i in range(nM):

            self.fobj[i] = Mode(nN)
            self.modalResults[i] = self.fobj[i].result




class Mode(object):


    def __init__(self, nN):

        self.freq = 0.0
        self.xmpf = 0.0
        self.ympf = 0.0
        self.zmpf = 0.0
        self.N = np.zeros(nN)
        self.xdsp = np.zeros(nN)
        self.ydsp = np.zeros(nN)
        self.zdsp = np.zeros(nN)
        self.xrot = np.zeros(nN)
        self.yrot = np.zeros(nN)
        self.zrot = np.zeros(nN)

        self.result = ModalResults(self.freq, self.xmpf, self.ympf, self.zmpf,
            ip(self.N), dp(self.xdsp), dp(self.ydsp), dp(self.zdsp),
            dp(self.xrot), dp(self.yrot), dp(self.zrot))






if __name__ == '__main__':


    N = np.arange(1, 13, dtype=np.int32)
    x = np.array([0.0, 120.0, 240.0, 360.0, 480.0, 600.0, 720.0, 120.0, 240.0, 360.0, 480.0, 600.0])
    y = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 120.0, 120.0, 120.0, 120.0, 120.0])
    z = np.zeros(12)
    r = np.zeros(12)

    # reactions
    NR = np.arange(1, 13, dtype=np.int32)
    Rx = np.array([1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0], dtype=np.int32)
    Ry = np.array([1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0], dtype=np.int32)
    Rz = np.ones(12, dtype=np.int32)
    Rxx = np.ones(12, dtype=np.int32)
    Ryy = np.ones(12, dtype=np.int32)
    Rzz = np.zeros(12, dtype=np.int32)

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

    # parameters
    shear = 0               # 1: include shear deformation
    geom = 0                # 1: include geometric stiffness
    exagg_static = 10.0     # exaggerate mesh deformations
    dx = 10.0               # x-axis increment for internal forces


    frame = Frame(N, x, y, z, r, NR, Rx, Ry, Rz, Rxx, Ryy, Rzz, EL, N1, N2, Ax, Asy, Asz,
        Jx, Iy, Iz, E, G, roll, density, shear, geom, exagg_static, dx)



    # load cases
    gx = 0.0
    gy = -386.4
    gz = 0.0

    load = StaticLoadCase(gx, gy, gz)

    nF = np.array([2, 3, 4, 5, 6], dtype=np.int32)
    Fx = np.zeros(5)
    Fy = np.array([-10.0, -20.0, -20.0, -10.0, -20.0])
    Fz = np.zeros(5)
    Mxx = np.zeros(5)
    Myy = np.zeros(5)
    Mzz = np.zeros(5)


    load.changePointLoads(nF, Fx, Fy, Fz, Mxx, Myy, Mzz)

    nD = np.array([8], dtype=np.int32)
    Dx = np.array([0.1])
    Dy = np.array([0.0])
    Dz = np.array([0.0])
    Dxx = np.array([0.0])
    Dyy = np.array([0.0])
    Dzz = np.array([0.0])

    load.changePrescribedDisplacements(nD, Dx, Dy, Dz, Dxx, Dyy, Dzz)


    frame.addLoadCase(load)
    frame.run()
