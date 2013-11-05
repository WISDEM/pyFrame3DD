#!/usr/bin/env python
# encoding: utf-8
"""
frame.py

Created by Andrew Ning on 2013-11-01.
Copyright (c) NREL. All rights reserved.
"""

import numpy as np
import math
from ctypes import POINTER, c_int, c_double, Structure, pointer
from collections import namedtuple


c_int_p = POINTER(c_int)
c_double_p = POINTER(c_double)


def ip(x):
    return x.ctypes.data_as(c_int_p)


def dp(x):
    return x.ctypes.data_as(c_double_p)



# --------------
# General Inputs
# --------------

class C_Nodes(Structure):
    _fields_ = [('nN', c_int),
                ('N', c_int_p),
                ('x', c_double_p),
                ('y', c_double_p),
                ('z', c_double_p),
                ('r', c_double_p)]



class C_Reactions(Structure):
    _fields_ = [('nR', c_int),
                ('N', c_int_p),
                ('Rx', c_int_p),
                ('Ry', c_int_p),
                ('Rz', c_int_p),
                ('Rxx', c_int_p),
                ('Ryy', c_int_p),
                ('Rzz', c_int_p)]



class C_Elements(Structure):
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



class C_OtherElementData(Structure):
    _fields_ = [('shear', c_int),
                ('geom', c_int),
                ('exagg_static', c_double),
                ('dx', c_double)]



# --------------
# Load Inputs
# --------------


class C_PointLoads(Structure):
    _fields_ = [('nF', c_int),
                ('N', c_int_p),
                ('Fx', c_double_p),
                ('Fy', c_double_p),
                ('Fz', c_double_p),
                ('Mxx', c_double_p),
                ('Myy', c_double_p),
                ('Mzz', c_double_p)]




class C_UniformLoads(Structure):
    _fields_ = [('nU', c_int),
                ('EL', c_int_p),
                ('Ux', c_double_p),
                ('Uy', c_double_p),
                ('Uz', c_double_p)]



class C_TrapezoidalLoads(Structure):
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



class C_ElementLoads(Structure):
    _fields_ = [('nP', c_int),
                ('EL', c_int_p),
                ('Px', c_double_p),
                ('Py', c_double_p),
                ('Pz', c_double_p),
                ('x', c_double_p)]


class C_TemperatureLoads(Structure):
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




class C_PrescribedDisplacements(Structure):
    _fields_ = [('nD', c_int),
                ('N', c_int_p),
                ('Dx', c_double_p),
                ('Dy', c_double_p),
                ('Dz', c_double_p),
                ('Dxx', c_double_p),
                ('Dyy', c_double_p),
                ('Dzz', c_double_p)]




class C_LoadCase(Structure):
    _fields_ = [('gx', c_double),
                ('gy', c_double),
                ('gz', c_double),
                ('pointLoads', C_PointLoads),
                ('uniformLoads', C_UniformLoads),
                ('trapezoidalLoads', C_TrapezoidalLoads),
                ('elementLoads', C_ElementLoads),
                ('temperatureLoads', C_TemperatureLoads),
                ('prescribedDisplacements', C_PrescribedDisplacements)]


# --------------
# Dynamic Inputs
# --------------


class C_DynamicData(Structure):
    _fields_ = [('nM', c_int),
                ('Mmethod', c_int),
                ('lump', c_int),
                ('tol', c_double),
                ('shift', c_double),
                ('exagg_modal', c_double)]


class C_ExtraInertia(Structure):
    _fields_ = [('nI', c_int),
                ('N', c_int_p),
                ('EMs', c_double_p),
                ('EMx', c_double_p),
                ('EMy', c_double_p),
                ('EMz', c_double_p)]


class C_ExtraMass(Structure):
    _fields_ = [('nX', c_int),
                ('EL', c_int_p),
                ('EMs', c_double_p)]


class C_Condensation(Structure):
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

class C_Displacements(Structure):
    _fields_ = [('node', c_int_p),
                ('x', c_double_p),
                ('y', c_double_p),
                ('z', c_double_p),
                ('xrot', c_double_p),
                ('yrot', c_double_p),
                ('zrot', c_double_p)]



class C_Forces(Structure):
    _fields_ = [('element', c_int_p),
                ('node', c_int_p),
                ('Nx', c_double_p),
                ('Vy', c_double_p),
                ('Vz', c_double_p),
                ('Txx', c_double_p),
                ('Myy', c_double_p),
                ('Mzz', c_double_p)]




class C_ReactionForces(Structure):
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


class C_InternalForces(Structure):
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


class C_MassResults(Structure):
    _fields_ = [('total_mass', c_double_p),
                ('struct_mass', c_double_p),
                ('N', c_int_p),
                ('xmass', c_double_p),
                ('ymass', c_double_p),
                ('zmass', c_double_p),
                ('xinrta', c_double_p),
                ('yinrta', c_double_p),
                ('zinrta', c_double_p),
                ]



class C_ModalResults(Structure):
    _fields_ = [('freq', c_double_p),
                ('xmpf', c_double_p),
                ('ympf', c_double_p),
                ('zmpf', c_double_p),
                ('N', c_int_p),
                ('xdsp', c_double_p),
                ('ydsp', c_double_p),
                ('zdsp', c_double_p),
                ('xrot', c_double_p),
                ('yrot', c_double_p),
                ('zrot', c_double_p),
                ]




# inputs

NodeData = namedtuple('NodeData', ['node', 'x', 'y', 'z', 'r'])
ReactionData = namedtuple('ReactionData', ['node', 'Rx', 'Ry', 'Rz', 'Rxx', 'Ryy', 'Rzz'])
ElementData = namedtuple('ElementData', ['element', 'N1', 'N2', 'Ax', 'Asy', 'Asz',
    'Jx', 'Iy', 'Iz', 'E', 'G', 'roll', 'density'])
OtherData = namedtuple('OtherData', ['shear', 'geom', 'exagg_static', 'dx'])


# outputs

NodeDisplacements = namedtuple('NodeDisplacements', ['node', 'dx', 'dy', 'dz', 'dxrot', 'dyrot', 'dzrot'])
ElementEndForces = namedtuple('ElementEndForces', ['element', 'node', 'Nx', 'Vy', 'Vz', 'Txx', 'Myy', 'Mzz'])
NodeReactions = namedtuple('NodeReactions', ['node', 'Fx', 'Fy', 'Fz', 'Mxx', 'Myy', 'Mzz'])
InternalForces = namedtuple('InternalForces', ['x', 'Nx', 'Vy', 'Vz', 'Tx', 'My', 'Mz', 'Dx', 'Dy', 'Dz', 'Rx'])
NodeMasses = namedtuple('NodeMasses', ['total_mass', 'struct_mass', 'node', 'xmass', 'ymass', 'zmass', 'xinrta', 'yinrta', 'zinrta'])
Modes = namedtuple('Modes', ['freq', 'xmpf', 'ympf', 'zmpf', 'node', 'xdsp', 'ydsp', 'zdsp',
    'xrot', 'yrot', 'zrot'])



class Frame(object):
    """docstring for Frame"""


    def __init__(self, nodes, reactions, elements, other):

        self.nodes = nodes
        self.reactions = reactions
        self.elements = elements
        self.other = other

        # convert to C int size (not longs) and copy to prevent releasing (b/c address space is shared by c)
        self.nnode = nodes.node.astype(np.int32)
        self.nx = np.copy(nodes.x)
        self.ny = np.copy(nodes.y)
        self.nz = np.copy(nodes.z)
        self.nr = np.copy(nodes.r)

        self.rnode = reactions.node.astype(np.int32)
        self.rRx = reactions.Rx.astype(np.int32)
        self.rRy = reactions.Ry.astype(np.int32)
        self.rRz = reactions.Rz.astype(np.int32)
        self.rRxx = reactions.Rxx.astype(np.int32)
        self.rRyy = reactions.Ryy.astype(np.int32)
        self.rRzz = reactions.Rzz.astype(np.int32)

        self.eelement = elements.element.astype(np.int32)
        self.eN1 = elements.N1.astype(np.int32)
        self.eN2 = elements.N2.astype(np.int32)
        self.eAx = np.copy(elements.Ax)
        self.eAsy = np.copy(elements.Asy)
        self.eAsz = np.copy(elements.Asz)
        self.eJx = np.copy(elements.Jx)
        self.eIy = np.copy(elements.Iy)
        self.eIz = np.copy(elements.Iz)
        self.eE = np.copy(elements.E)
        self.eG = np.copy(elements.G)
        self.eroll = np.copy(elements.roll)
        self.edensity = np.copy(elements.density)


        # create c objects
        self.c_nodes = C_Nodes(len(self.nnode), ip(self.nnode), dp(self.nx),
            dp(self.ny), dp(self.nz), dp(self.nr))

        self.c_reactions = C_Reactions(len(self.rnode), ip(self.rnode),
            ip(self.rRx), ip(self.rRy), ip(self.rRz),
            ip(self.rRxx), ip(self.rRyy), ip(self.rRzz))

        self.c_elements = C_Elements(len(self.eelement), ip(self.eelement),
            ip(self.eN1), ip(self.eN2), dp(self.eAx), dp(self.eAsy),
            dp(self.eAsz), dp(self.eJx), dp(self.eIy), dp(self.eIz),
            dp(self.eE), dp(self.eG), dp(self.eroll), dp(self.edensity))

        self.c_other = C_OtherElementData(other.shear, other.geom, other.exagg_static, other.dx)


        # create list for load cases
        self.loadCases = []


        # initialize with no dynamics
        dynamic = DynamicAnalysis(nM=0, Mmethod=1, lump=0, tol=0.0, shift=0.0, exagg_modal=0.0)
        self.useDynamicAnalysis(dynamic)


        # load c module
        self._frame3dd = np.ctypeslib.load_library('_pyframe3dd', '.')

        self._frame3dd.run.argtypes = [POINTER(C_Nodes), POINTER(C_Reactions), POINTER(C_Elements),
            POINTER(C_OtherElementData), c_int, POINTER(C_LoadCase),
            POINTER(C_DynamicData), POINTER(C_ExtraInertia), POINTER(C_ExtraMass),
            POINTER(C_Condensation),
            POINTER(C_Displacements), POINTER(C_Forces), POINTER(C_ReactionForces),
            POINTER(POINTER(C_InternalForces)), POINTER(C_MassResults), POINTER(C_ModalResults)]

        self._frame3dd.run.restype = c_int



    def addLoadCase(self, loadCase):

        self.loadCases.append(loadCase)



    def useDynamicAnalysis(self, dynamic):

        self.dynamic = dynamic


    def run(self):

        nCases = len(self.loadCases)  # number of load cases
        nN = len(self.nodes.node)  # number of nodes
        nE = len(self.elements.element)  # number of elements
        nM = self.dynamic.nM  # number of modes

        if nCases == 0:
            print('error: must have at least 1 load case')
            return



        # initialize output arrays

        dout = NodeDisplacements(np.zeros((nCases, nN), dtype=np.int32),
            np.zeros((nCases, nN)), np.zeros((nCases, nN)), np.zeros((nCases, nN)),
            np.zeros((nCases, nN)), np.zeros((nCases, nN)), np.zeros((nCases, nN))
        )
        fout = ElementEndForces(np.zeros((nCases, 2*nE), dtype=np.int32),
            np.zeros((nCases, 2*nE), dtype=np.int32),
            np.zeros((nCases, 2*nE)), np.zeros((nCases, 2*nE)), np.zeros((nCases, 2*nE)),
            np.zeros((nCases, 2*nE)), np.zeros((nCases, 2*nE)), np.zeros((nCases, 2*nE))
        )
        rout = NodeReactions(np.zeros((nCases, nN), dtype=np.int32),
            np.zeros((nCases, nN)), np.zeros((nCases, nN)), np.zeros((nCases, nN)),
            np.zeros((nCases, nN)), np.zeros((nCases, nN)), np.zeros((nCases, nN))
        )

        x = self.nodes.x
        y = self.nodes.y
        z = self.nodes.z
        N1 = self.elements.N1
        N2 = self.elements.N2
        dx = self.other.dx

        ifout = [0]*nE
        for i in range(nE):

            L = math.sqrt(
                (x[N2[i]-1] - x[N1[i]-1])**2 +
                (y[N2[i]-1] - y[N1[i]-1])**2 +
                (z[N2[i]-1] - z[N1[i]-1])**2
            )

            nIF = int(max(math.floor(L/dx), 1)) + 1

            ifout[i] = InternalForces(np.zeros((nCases, nIF)), np.zeros((nCases, nIF)),
                np.zeros((nCases, nIF)), np.zeros((nCases, nIF)), np.zeros((nCases, nIF)),
                np.zeros((nCases, nIF)), np.zeros((nCases, nIF)), np.zeros((nCases, nIF)),
                np.zeros((nCases, nIF)), np.zeros((nCases, nIF)), np.zeros((nCases, nIF))
            )


        mout = NodeMasses(0.0, 0.0, np.zeros(nN, dtype=np.int32),
            np.zeros(nN), np.zeros(nN), np.zeros(nN),
            np.zeros(nN), np.zeros(nN), np.zeros(nN)
        )
        modalout = Modes(np.zeros(nM), np.zeros(nM), np.zeros(nM), np.zeros(nM),
            np.zeros((nM, nN), dtype=np.int32), np.zeros((nM, nN)), np.zeros((nM, nN)),
            np.zeros((nM, nN)), np.zeros((nM, nN)), np.zeros((nM, nN)),
            np.zeros((nM, nN))
        )



        # create c structs

        c_loadcases = (C_LoadCase * nCases)()
        c_disp = (C_Displacements * nCases)()
        c_forces = (C_Forces * nCases)()
        c_reactions = (C_ReactionForces * nCases)()
        c_internalForces = (POINTER(C_InternalForces) * nCases)()


        for i in range(nCases):
            lci = self.loadCases[i]
            c_loadcases[i] = C_LoadCase(lci.gx, lci.gy, lci.gz, lci.pL,
                lci.uL, lci.tL, lci.eL, lci.tempL, lci.pD)
            c_disp[i] = C_Displacements(ip(dout.node[i, :]),
                dp(dout.dx[i, :]), dp(dout.dy[i, :]), dp(dout.dz[i, :]),
                dp(dout.dxrot[i, :]), dp(dout.dyrot[i, :]), dp(dout.dzrot[i, :]))
            c_forces[i] = C_Forces(ip(fout.element[i, :]), ip(fout.node[i, :]),
                dp(fout.Nx[i, :]), dp(fout.Vy[i, :]), dp(fout.Vz[i, :]),
                dp(fout.Txx[i, :]), dp(fout.Myy[i, :]), dp(fout.Mzz[i, :]))
            c_reactions[i] = C_ReactionForces(ip(rout.node[i, :]),
                dp(rout.Fx[i, :]), dp(rout.Fy[i, :]), dp(rout.Fz[i, :]),
                dp(rout.Mxx[i, :]), dp(rout.Myy[i, :]), dp(rout.Mzz[i, :]))

            c_internalForces[i] = (C_InternalForces * nE)()
            for j in range(nE):
                (c_internalForces[i])[j] = C_InternalForces(dp(ifout[j].x[i, :]), dp(ifout[j].Nx[i, :]),
                    dp(ifout[j].Vy[i, :]), dp(ifout[j].Vz[i, :]), dp(ifout[j].Tx[i, :]),
                    dp(ifout[j].My[i, :]), dp(ifout[j].Mz[i, :]), dp(ifout[j].Dx[i, :]),
                    dp(ifout[j].Dy[i, :]), dp(ifout[j].Dz[i, :]), dp(ifout[j].Rx[i, :]))

        total_mass = c_double()
        struct_mass = c_double()

        c_massResults = C_MassResults(pointer(total_mass), pointer(struct_mass), ip(mout.node),
            dp(mout.xmass), dp(mout.ymass), dp(mout.zmass),
            dp(mout.xinrta), dp(mout.yinrta), dp(mout.zinrta))

        c_modalResults = (C_ModalResults * nM)()

        freq = [0]*nM
        xmpf = [0]*nM
        ympf = [0]*nM
        zmpf = [0]*nM

        for i in range(nM):

            freq[i] = c_double()
            xmpf[i] = c_double()
            ympf[i] = c_double()
            zmpf[i] = c_double()

            c_modalResults[i] = C_ModalResults(pointer(freq[i]), pointer(xmpf[i]),
                pointer(ympf[i]), pointer(zmpf[i]), ip(modalout.node[i, :]),
                dp(modalout.xdsp[i, :]), dp(modalout.ydsp[i, :]), dp(modalout.zdsp[i, :]),
                dp(modalout.xrot[i, :]), dp(modalout.yrot[i, :]), dp(modalout.zrot[i, :])
            )


        d = self.dynamic


        self._frame3dd.run(self.c_nodes, self.c_reactions, self.c_elements, self.c_other,
            nCases, c_loadcases,
            d.dynamicData, d.extraInertia, d.extraMass, d.condensation,
            c_disp, c_forces, c_reactions, c_internalForces, c_massResults, c_modalResults)

        # put mass values back in since tuple is read only
        mout = NodeMasses(total_mass.value, struct_mass.value, mout.node,
            mout.xmass, mout.ymass, mout.zmass,
            mout.xinrta, mout.yinrta, mout.zinrta)

        # put modal results back in
        for i in range(nM):
            modalout.freq[i] = freq[i].value
            modalout.xmpf[i] = xmpf[i].value
            modalout.ympf[i] = ympf[i].value
            modalout.zmpf[i] = zmpf[i].value

        return dout, fout, rout, ifout, mout, modalout





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

        # copying to prevent any user error with variables pointing to something else (b/c memory address is shared by C)
        self.NF = N.astype(np.int32)
        self.Fx = np.copy(Fx)
        self.Fy = np.copy(Fy)
        self.Fz = np.copy(Fz)
        self.Mxx = np.copy(Mxx)
        self.Myy = np.copy(Myy)
        self.Mzz = np.copy(Mzz)

        self.pL = C_PointLoads(len(N), ip(self.NF), dp(self.Fx), dp(self.Fy), dp(self.Fz),
                   dp(self.Mxx), dp(self.Myy), dp(self.Mzz))


    def changeUniformLoads(self, EL, Ux, Uy, Uz):

        self.ELU = EL.astype(np.int32)
        self.Ux = np.copy(Ux)
        self.Uy = np.copy(Uy)
        self.Uz = np.copy(Uz)

        self.uL = C_UniformLoads(len(EL), ip(self.ELU), dp(self.Ux), dp(self.Uy), dp(self.Uz))



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

        self.tL = C_TrapezoidalLoads(len(EL), ip(self.ELT), dp(self.xx1), dp(self.xx2), dp(self.wx1), dp(self.wx2),
            dp(self.xy1), dp(self.xy2), dp(self.wy1), dp(self.wy2), dp(self.xz1), dp(self.xz2), dp(self.wz1), dp(self.wz2))



    def changeElementLoads(self, EL, Px, Py, Pz, x):

        self.ELE = EL.astype(np.int32)
        self.Px = np.copy(Px)
        self.Py = np.copy(Py)
        self.Pz = np.copy(Pz)
        self.xE = np.copy(x)


        self.eL = C_ElementLoads(len(EL), ip(self.ELE), dp(self.Px), dp(self.Py), dp(self.Pz), dp(self.xE))



    def changeTemperatureLoads(self, EL, a, hy, hz, Typ, Tym, Tzp, Tzm):

        self.ELTemp = EL.astype(np.int32)
        self.a = np.copy(a)
        self.hy = np.copy(hy)
        self.hz = np.copy(hz)
        self.Typ = np.copy(Typ)
        self.Tym = np.copy(Tym)
        self.Tzp = np.copy(Tzp)
        self.Tzm = np.copy(Tzm)


        self.tempL = C_TemperatureLoads(len(EL), ip(self.ELTemp), dp(self.a), dp(self.hy), dp(self.hz), dp(self.Typ),
            dp(self.Tym), dp(self.Tzp), dp(self.Tzm))



    def changePrescribedDisplacements(self, N, Dx, Dy, Dz, Dxx, Dyy, Dzz):

        self.ND = N.astype(np.int32)
        self.Dx = np.copy(Dx)
        self.Dy = np.copy(Dy)
        self.Dz = np.copy(Dz)
        self.Dxx = np.copy(Dxx)
        self.Dyy = np.copy(Dyy)
        self.Dzz = np.copy(Dzz)

        self.pD = C_PrescribedDisplacements(len(N), ip(self.ND), dp(self.Dx), dp(self.Dy), dp(self.Dz),
                   dp(self.Dxx), dp(self.Dyy), dp(self.Dzz))




class DynamicAnalysis(object):
    """docstring"""


    def __init__(self, nM, Mmethod, lump, tol, shift, exagg_modal):

        self.nM = nM

        self.dynamicData = C_DynamicData(nM, Mmethod, lump, tol, shift, exagg_modal)

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


        self.extraInertia = C_ExtraInertia(len(N), ip(self.NI), dp(self.EMs),
            dp(self.EMx), dp(self.EMy), dp(self.EMz))



    def changeExtraMass(self, EL, EMs):

        self.ELM = EL.astype(np.int32)
        self.EMsM = np.copy(EMs)

        self.extraMass = C_ExtraMass(len(EL), ip(self.ELM), dp(self.EMsM))




    def changeCondensationData(self, Cmethod, N, cx, cy, cz, cxx, cyy, czz, m):

        self.NC = N.astype(np.int32)
        self.cx = np.copy(cx)
        self.cy = np.copy(cy)
        self.cz = np.copy(cx)
        self.cxx = np.copy(cxx)
        self.cyy = np.copy(cyy)
        self.czz = np.copy(czz)
        self.mC = m.astype(np.int32)

        self.condensation = C_Condensation(Cmethod, len(N), ip(self.NC), dp(self.cx), dp(self.cy), dp(self.cz),
            dp(self.cxx), dp(self.cyy), dp(self.czz), ip(self.mC))







if __name__ == '__main__':


    # nodes
    node = np.arange(1, 13)
    x = np.array([0.0, 120.0, 240.0, 360.0, 480.0, 600.0, 720.0, 120.0, 240.0, 360.0, 480.0, 600.0])
    y = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 120.0, 120.0, 120.0, 120.0, 120.0])
    z = np.zeros(12)
    r = np.zeros(12)
    nodes = NodeData(node, x, y, z, r)

    # reactions
    node = np.arange(1, 13)
    Rx = np.array([1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0])
    Ry = np.array([1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0])
    Rz = np.ones(12)
    Rxx = np.ones(12)
    Ryy = np.ones(12)
    Rzz = np.zeros(12)
    reactions = ReactionData(node, Rx, Ry, Rz, Rxx, Ryy, Rzz)

    # elements
    EL = np.arange(1, 22)
    N1 = np.array([1, 2, 3, 4, 5, 6, 1, 2, 2, 3, 4, 4, 4, 5, 6, 6, 7, 8, 9, 10, 11])
    N2 = np.array([2, 3, 4, 5, 6, 7, 8, 8, 9, 9, 9, 10, 11, 11, 11, 12, 12, 9, 10, 11, 12])
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
    elements = ElementData(EL, N1, N2, Ax, Asy, Asz, Jx, Iy, Iz, E, G, roll, density)

    # parameters
    shear = 0               # 1: include shear deformation
    geom = 0                # 1: include geometric stiffness
    exagg_static = 10.0     # exaggerate mesh deformations
    dx = 10.0               # x-axis increment for internal forces
    other = OtherData(shear, geom, exagg_static, dx)


    frame = Frame(nodes, reactions, elements, other)



    # load cases
    gx = 0.0
    gy = -386.4
    gz = 0.0

    load = StaticLoadCase(gx, gy, gz)

    nF = np.array([2, 3, 4, 5, 6])
    Fx = np.zeros(5)
    Fy = np.array([-10.0, -20.0, -20.0, -10.0, -20.0])
    Fz = np.zeros(5)
    Mxx = np.zeros(5)
    Myy = np.zeros(5)
    Mzz = np.zeros(5)


    load.changePointLoads(nF, Fx, Fy, Fz, Mxx, Myy, Mzz)

    nD = np.array([8])
    Dx = np.array([0.1])
    Dy = np.array([0.0])
    Dz = np.array([0.0])
    Dxx = np.array([0.0])
    Dyy = np.array([0.0])
    Dzz = np.array([0.0])

    load.changePrescribedDisplacements(nD, Dx, Dy, Dz, Dxx, Dyy, Dzz)


    frame.addLoadCase(load)
    displacements, forces, reactions, internalForces, mass, modal = frame.run()

    iCase = 0

    print 'nodes =', displacements.node[iCase, :]
    print 'dx =', displacements.dx[iCase, :]
    print 'dy =', displacements.dy[iCase, :]
    print 'dz =', displacements.dz[iCase, :]
    print 'dxrot =', displacements.dxrot[iCase, :]
    print 'dyrot =', displacements.dyrot[iCase, :]
    print 'dzrot =', displacements.dzrot[iCase, :]

    print

    print 'element =', forces.element[iCase, :]
    print 'node =', forces.node[iCase, :]
    print 'Nx =', forces.Nx[iCase, :]
    print 'Vy =', forces.Vy[iCase, :]
    print 'Vz =', forces.Vz[iCase, :]
    print 'Txx =', forces.Txx[iCase, :]
    print 'Myy =', forces.Myy[iCase, :]
    print 'Mzz =', forces.Mzz[iCase, :]

    print

    print 'nodesR =', reactions.node[iCase, :]
    print 'RFx =', reactions.Fx[iCase, :]
    print 'RFy =', reactions.Fy[iCase, :]
    print 'RFz =', reactions.Fz[iCase, :]
    print 'RMxx =', reactions.Mxx[iCase, :]
    print 'RMyy =', reactions.Myy[iCase, :]
    print 'RMzz =', reactions.Mzz[iCase, :]

    print


    iE = 2

    print 'xE =', internalForces[iE].x[iCase, :]
    print 'Nx =', internalForces[iE].Nx[iCase, :]
    print 'Vy =', internalForces[iE].Vy[iCase, :]
    print 'Vz =', internalForces[iE].Vz[iCase, :]
    print 'Tx =', internalForces[iE].Tx[iCase, :]
    print 'My =', internalForces[iE].My[iCase, :]
    print 'Mz =', internalForces[iE].Mz[iCase, :]
    print 'Dx =', internalForces[iE].Dx[iCase, :]
    print 'Dy =', internalForces[iE].Dy[iCase, :]
    print 'Dz =', internalForces[iE].Dz[iCase, :]
    print 'Rx =', internalForces[iE].Rx[iCase, :]

    print

    print mass.total_mass
    print mass.struct_mass

    print

    print modal.freq

