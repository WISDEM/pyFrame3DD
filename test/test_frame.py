#!/usr/bin/env python
# encoding: utf-8
"""
untitled.py

Created by Andrew Ning on 2013-11-04.
Copyright (c) NREL. All rights reserved.
"""

import unittest
import numpy as np
from StringIO import StringIO

from frame3dd import Frame, NodeData, ReactionData, ElementData, OtherData, \
    StaticLoadCase, DynamicAnalysis


class FrameTestEXA(unittest.TestCase):

    def setUp(self):

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



        # load cases 1
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


        # load cases
        gx = 0.0
        gy = -386.4
        gz = 0.0

        load = StaticLoadCase(gx, gy, gz)

        nF = np.array([3, 4, 5])
        Fx = np.array([20.0, 10.0, 20.0])
        Fy = np.zeros(3)
        Fz = np.zeros(3)
        Mxx = np.zeros(3)
        Myy = np.zeros(3)
        Mzz = np.zeros(3)

        load.changePointLoads(nF, Fx, Fy, Fz, Mxx, Myy, Mzz)

        EL = np.array([10, 13, 15])
        a = 6e-12 * np.ones(3)
        hy = 5.0 * np.ones(3)
        hz = 5.0 * np.ones(3)
        Typ = np.array([10.0, 15.0, 17.0])
        Tym = np.array([10.0, 15.0, 17.0])
        Tzp = np.array([10.0, 15.0, 17.0])
        Tzm = np.array([10.0, 15.0, 17.0])

        load.changeTemperatureLoads(EL, a, hy, hz, Typ, Tym, Tzp, Tzm)


        nD = np.array([1, 8])
        Dx = np.array([0.0, 0.1])
        Dy = np.array([-1.0, 0.0])
        Dz = np.array([0.0, 0.0])
        Dxx = np.array([0.0, 0.0])
        Dyy = np.array([0.0, 0.0])
        Dzz = np.array([0.0, 0.0])

        load.changePrescribedDisplacements(nD, Dx, Dy, Dz, Dxx, Dyy, Dzz)


        frame.addLoadCase(load)

        self.displacements, self.forces, self.reactions, self.internalForces, self.mass, self.modal = frame.run()






    def test_disp1(self):

        disp = self.displacements
        iCase = 0

        node = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], dtype=np.int)
        dx = np.array([0.0, 0.010776, 0.035528, 0.060279, 0.086295, 0.112311, 0.129754, 0.100000, 0.089226, 0.059394, 0.029563, 0.012122])
        dy = np.array([0.0, -0.171344, -0.297816, -0.332643, -0.295487, -0.184135, 0.0, -0.152896, -0.289325, -0.332855, -0.291133, -0.166951])
        dz = np.zeros(12)
        dxrot = np.zeros(12)
        dyrot = np.zeros(12)
        dzrot = np.array([-0.501728, -0.087383, 0.015516, 0.000017, -0.015556, 0.087387, 0.501896, 0.136060, -0.011102, -0.000000, 0.011062, -0.136029])

        np.testing.assert_array_equal(disp.node[iCase, :], node)
        np.testing.assert_array_almost_equal(disp.dx[iCase, :], dx, decimal=6)
        np.testing.assert_array_almost_equal(disp.dy[iCase, :], dy, decimal=6)
        np.testing.assert_array_almost_equal(disp.dz[iCase, :], dz, decimal=6)
        np.testing.assert_array_almost_equal(disp.dxrot[iCase, :], dxrot, decimal=6)
        np.testing.assert_array_almost_equal(disp.dyrot[iCase, :], dyrot, decimal=6)
        np.testing.assert_array_almost_equal(disp.dzrot[iCase, :], dzrot, decimal=6)


    def test_disp2(self):

        disp = self.displacements
        iCase = 1

        node = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], dtype=np.int)
        dx = np.array([0.0, 0.071965, 0.134909, 0.189577, 0.220207, 0.242560, 0.254035, 0.100000, 0.048727, 0.016113, -0.016502, -0.027973])
        dy = np.array([-1.000000, -1.067463, -1.018927, -0.850595, -0.615710, -0.325659, 0.0, -1.076148, -1.018711, -0.850807, -0.615495, -0.314444])
        dz = np.zeros(12)
        dxrot = np.zeros(12)
        dyrot = np.zeros(12)
        dzrot = np.array([-0.501206, -0.086438, 0.016991, 0.001616, -0.013988, 0.088765, 0.503041, 0.136834, -0.009551, 0.001627, 0.012588, -0.134603])

        np.testing.assert_array_equal(disp.node[iCase, :], node)
        np.testing.assert_array_almost_equal(disp.dx[iCase, :], dx, decimal=6)
        np.testing.assert_array_almost_equal(disp.dy[iCase, :], dy, decimal=6)
        np.testing.assert_array_almost_equal(disp.dz[iCase, :], dz, decimal=6)
        np.testing.assert_array_almost_equal(disp.dxrot[iCase, :], dxrot, decimal=6)
        np.testing.assert_array_almost_equal(disp.dyrot[iCase, :], dyrot, decimal=6)
        np.testing.assert_array_almost_equal(disp.dzrot[iCase, :], dzrot, decimal=6)


    def test_force1(self):

        forces = self.forces
        iCase = 0

        output = StringIO("""
         1      1    -26.042      0.099      0.0        0.0        0.0       -1.853
         1      2     26.042      0.241      0.0        0.0        0.0       -6.648
         2      2    -59.816      0.162      0.0        0.0        0.0        2.644
         2      3     59.816      0.178      0.0        0.0        0.0       -3.656
         3      3    -59.815      0.172      0.0        0.0        0.0        3.553
         3      4     59.815      0.168      0.0        0.0        0.0       -3.319
         4      4    -62.872      0.168      0.0        0.0        0.0        3.319
         4      5     62.872      0.172      0.0        0.0        0.0       -3.554
         5      5    -62.872      0.178      0.0        0.0        0.0        3.657
         5      6     62.872      0.161      0.0        0.0        0.0       -2.643
         6      6    -42.155      0.241      0.0        0.0        0.0        6.647
         6      7     42.155      0.099      0.0        0.0        0.0        1.853
         7      1     64.086      0.148      0.0        0.0        0.0        1.853
         7      8    -63.746      0.192      0.0        0.0        0.0       -5.581
         8      2    -44.414      0.006      0.0        0.0        0.0       -0.176
         8      8     44.754     -0.006      0.0        0.0        0.0        0.904
         9      2     47.936      0.164      0.0        0.0        0.0        4.180
         9      9    -47.596      0.176      0.0        0.0        0.0       -5.173
        10      3    -20.350      0.001      0.0        0.0        0.0        0.103
        10      9     20.690     -0.001      0.0        0.0        0.0       -0.026
        11      4    -17.194     -0.171      0.0        0.0        0.0       -4.841
        11      9     17.534     -0.169      0.0        0.0        0.0        4.734
        12      4      0.682      0.0        0.0        0.0        0.0        0.0
        12     10     -0.342      0.0        0.0        0.0        0.0        0.0
        13      4    -12.872      0.171      0.0        0.0        0.0        4.841
        13     11     13.212      0.169      0.0        0.0        0.0       -4.734
        14      5    -10.350     -0.001      0.0        0.0        0.0       -0.104
        14     11     10.690      0.001      0.0        0.0        0.0        0.025
        15      6     29.472     -0.164      0.0        0.0        0.0       -4.180
        15     11    -29.132     -0.176      0.0        0.0        0.0        5.173
        16      6    -41.358     -0.006      0.0        0.0        0.0        0.175
        16     12     41.698      0.006      0.0        0.0        0.0       -0.905
        17      7     59.764     -0.148      0.0        0.0        0.0       -1.853
        17     12    -59.424     -0.192      0.0        0.0        0.0        5.580
        18      8     26.036      0.185      0.0        0.0        0.0        4.677
        18      9    -26.036      0.155      0.0        0.0        0.0       -2.832
        19      9     72.094      0.169      0.0        0.0        0.0        3.297
        19     10    -72.094      0.171      0.0        0.0        0.0       -3.447
        20     10     72.094      0.171      0.0        0.0        0.0        3.447
        20     11    -72.094      0.169      0.0        0.0        0.0       -3.297
        21     11     42.149      0.155      0.0        0.0        0.0        2.833
        21     12    -42.149      0.185      0.0        0.0        0.0       -4.675
        """)

        out = np.loadtxt(output)

        np.testing.assert_array_equal(forces.element[iCase, :], out[:, 0])
        np.testing.assert_array_equal(forces.node[iCase, :], out[:, 1])
        np.testing.assert_array_almost_equal(forces.Nx[iCase, :], out[:, 2], decimal=3)
        np.testing.assert_array_almost_equal(forces.Vy[iCase, :], out[:, 3], decimal=3)
        np.testing.assert_array_almost_equal(forces.Vz[iCase, :], out[:, 4], decimal=3)
        np.testing.assert_array_almost_equal(forces.Txx[iCase, :], out[:, 5], decimal=3)
        np.testing.assert_array_almost_equal(forces.Myy[iCase, :], out[:, 6], decimal=3)
        np.testing.assert_array_almost_equal(forces.Mzz[iCase, :], out[:, 7], decimal=3)


    def test_force2(self):

        forces = self.forces
        iCase = 1

        output = StringIO("""
         1      1   -173.916      0.099      0.0        0.0        0.0       -1.856
         1      2    173.916      0.241      0.0        0.0        0.0       -6.649
         2      2   -152.115      0.161      0.0        0.0        0.0        2.639
         2      3    152.115      0.178      0.0        0.0        0.0       -3.658
         3      3   -132.114      0.172      0.0        0.0        0.0        3.550
         3      4    132.114      0.168      0.0        0.0        0.0       -3.321
         4      4    -74.021      0.168      0.0        0.0        0.0        3.318
         4      5     74.021      0.172      0.0        0.0        0.0       -3.555
         5      5    -54.022      0.178      0.0        0.0        0.0        3.658
         5      6     54.022      0.161      0.0        0.0        0.0       -2.643
         6      6    -27.729      0.241      0.0        0.0        0.0        6.649
         6      7     27.729      0.099      0.0        0.0        0.0        1.854
         7      1    -28.651      0.148      0.0        0.0        0.0        1.856
         7      8     28.991      0.192      0.0        0.0        0.0       -5.577
         8      2     21.160      0.006      0.0        0.0        0.0       -0.171
         8      8    -20.820     -0.006      0.0        0.0        0.0        0.908
         9      2    -30.658      0.164      0.0        0.0        0.0        4.180
         9      9     30.998      0.176      0.0        0.0        0.0       -5.170
        10      3     -0.350      0.001      0.0        0.0        0.0        0.108
        10      9      0.690     -0.001      0.0        0.0        0.0       -0.021
        11      4     33.116     -0.171      0.0        0.0        0.0       -4.841
        11      9    -32.776     -0.169      0.0        0.0        0.0        4.734
        12      4      0.682      0.0        0.0        0.0        0.0        0.003
        12     10     -0.342      0.0        0.0        0.0        0.0        0.003
        13      4    -34.898      0.171      0.0        0.0        0.0        4.842
        13     11     35.237      0.169      0.0        0.0        0.0       -4.734
        14      5     -0.350     -0.001      0.0        0.0        0.0       -0.103
        14     11      0.690      0.001      0.0        0.0        0.0        0.025
        15      6     37.356     -0.164      0.0        0.0        0.0       -4.180
        15     11    -37.016     -0.176      0.0        0.0        0.0        5.173
        16      6    -26.933     -0.006      0.0        0.0        0.0        0.175
        16     12     27.273      0.006      0.0        0.0        0.0       -0.905
        17      7     39.363     -0.148      0.0        0.0        0.0       -1.854
        17     12    -39.023     -0.192      0.0        0.0        0.0        5.580
        18      8    123.910      0.185      0.0        0.0        0.0        4.668
        18      9   -123.910      0.155      0.0        0.0        0.0       -2.837
        19      9     78.818      0.169      0.0        0.0        0.0        3.294
        19     10    -78.818      0.171      0.0        0.0        0.0       -3.449
        20     10     78.818      0.171      0.0        0.0        0.0        3.447
        20     11    -78.818      0.169      0.0        0.0        0.0       -3.298
        21     11     27.723      0.155      0.0        0.0        0.0        2.834
        21     12    -27.723      0.185      0.0        0.0        0.0       -4.675
        """)

        out = np.loadtxt(output)

        np.testing.assert_array_equal(forces.element[iCase, :], out[:, 0])
        np.testing.assert_array_equal(forces.node[iCase, :], out[:, 1])
        np.testing.assert_array_almost_equal(forces.Nx[iCase, :], out[:, 2], decimal=3)
        np.testing.assert_array_almost_equal(forces.Vy[iCase, :], out[:, 3], decimal=3)
        np.testing.assert_array_almost_equal(forces.Vz[iCase, :], out[:, 4], decimal=3)
        np.testing.assert_array_almost_equal(forces.Txx[iCase, :], out[:, 5], decimal=3)
        np.testing.assert_array_almost_equal(forces.Myy[iCase, :], out[:, 6], decimal=3)
        np.testing.assert_array_almost_equal(forces.Mzz[iCase, :], out[:, 7], decimal=3)



    def test_reactions1(self):

        reactions = self.reactions
        iCase = 0

        output = StringIO("""
             1      19.168      45.519       0.0         0.0         0.0         0.0
             2       0.0         0.0         0.0         0.0         0.0         0.0
             3       0.0         0.0         0.0         0.0         0.0         0.0
             4       0.0         0.0         0.0         0.0         0.0         0.0
             5       0.0         0.0         0.0         0.0         0.0         0.0
             6       0.0         0.0         0.0         0.0         0.0         0.0
             7       0.0        42.463       0.0         0.0         0.0         0.0
             8     -19.168       0.0         0.0         0.0         0.0         0.0
             9       0.0         0.0         0.0         0.0         0.0         0.0
            10       0.0         0.0         0.0         0.0         0.0         0.0
            11       0.0         0.0         0.0         0.0         0.0         0.0
            12       0.0         0.0         0.0         0.0         0.0         0.0
        """)

        out = np.loadtxt(output)

        np.testing.assert_array_equal(reactions.node[iCase, :], out[:, 0])
        np.testing.assert_array_almost_equal(reactions.Fx[iCase, :], out[:, 1], decimal=3)
        np.testing.assert_array_almost_equal(reactions.Fy[iCase, :], out[:, 2], decimal=3)
        np.testing.assert_array_almost_equal(reactions.Fz[iCase, :], out[:, 3], decimal=3)
        np.testing.assert_array_almost_equal(reactions.Mxx[iCase, :], out[:, 4], decimal=3)
        np.testing.assert_array_almost_equal(reactions.Myy[iCase, :], out[:, 5], decimal=3)
        np.testing.assert_array_almost_equal(reactions.Mzz[iCase, :], out[:, 6], decimal=3)





    def test_reactions2(self):

        reactions = self.reactions
        iCase = 1

        output = StringIO("""
         1    -194.280     -20.056       0.0         0.0         0.0         0.0
         2       0.0         0.0         0.0         0.0         0.0         0.0
         3       0.0         0.0         0.0         0.0         0.0         0.0
         4       0.0         0.0         0.0         0.0         0.0         0.0
         5       0.0         0.0         0.0         0.0         0.0         0.0
         6       0.0         0.0         0.0         0.0         0.0         0.0
         7       0.0        28.038       0.0         0.0         0.0         0.0
         8     144.280       0.0         0.0         0.0         0.0         0.0
         9       0.0         0.0         0.0         0.0         0.0         0.0
        10       0.0         0.0         0.0         0.0         0.0         0.0
        11       0.0         0.0         0.0         0.0         0.0         0.0
        12       0.0         0.0         0.0         0.0         0.0         0.0
        """)

        out = np.loadtxt(output)

        np.testing.assert_array_equal(reactions.node[iCase, :], out[:, 0])
        np.testing.assert_array_almost_equal(reactions.Fx[iCase, :], out[:, 1], decimal=3)
        np.testing.assert_array_almost_equal(reactions.Fy[iCase, :], out[:, 2], decimal=3)
        np.testing.assert_array_almost_equal(reactions.Fz[iCase, :], out[:, 3], decimal=3)
        np.testing.assert_array_almost_equal(reactions.Mxx[iCase, :], out[:, 4], decimal=3)
        np.testing.assert_array_almost_equal(reactions.Myy[iCase, :], out[:, 5], decimal=3)
        np.testing.assert_array_almost_equal(reactions.Mzz[iCase, :], out[:, 6], decimal=3)



    def test_if1(self):

        intF = self.internalForces
        iE = 3
        iCase = 0

        output = StringIO("""
            0.000000e+00    6.287160e+01   -1.679863e-01    0.000000e+00    0.000000e+00    0.000000e+00   -3.319261e+00    6.027882e-02   -3.326427e-01    0.000000e+00    0.000000e+00
            1.000000e+01    6.287160e+01   -1.396631e-01    0.000000e+00    0.000000e+00    0.000000e+00   -1.781014e+00    6.244680e-02   -7.675230e-01    0.000000e+00    0.000000e+00
            2.000000e+01    6.287160e+01   -1.113400e-01    0.000000e+00    0.000000e+00    0.000000e+00   -5.259979e-01    6.461479e-02   -1.832824e+00    0.000000e+00    0.000000e+00
            3.000000e+01    6.287160e+01   -8.301690e-02    0.000000e+00    0.000000e+00    0.000000e+00    4.457867e-01    6.678278e-02   -3.095780e+00    0.000000e+00    0.000000e+00
            4.000000e+01    6.287160e+01   -5.469378e-02    0.000000e+00    0.000000e+00    0.000000e+00    1.134340e+00    6.895076e-02   -4.221295e+00    0.000000e+00    0.000000e+00
            5.000000e+01    6.287160e+01   -2.637066e-02    0.000000e+00    0.000000e+00    0.000000e+00    1.539662e+00    7.111875e-02   -4.971936e+00    0.000000e+00    0.000000e+00
            6.000000e+01    6.287160e+01    1.952455e-03    0.000000e+00    0.000000e+00    0.000000e+00    1.661753e+00    7.328673e-02   -5.207937e+00    0.000000e+00    0.000000e+00
            7.000000e+01    6.287160e+01    3.027557e-02    0.000000e+00    0.000000e+00    0.000000e+00    1.500613e+00    7.545472e-02   -4.887197e+00    0.000000e+00    0.000000e+00
            8.000000e+01    6.287160e+01    5.859869e-02    0.000000e+00    0.000000e+00    0.000000e+00    1.056242e+00    7.762271e-02   -4.065281e+00    0.000000e+00    0.000000e+00
            9.000000e+01    6.287160e+01    8.692181e-02    0.000000e+00    0.000000e+00    0.000000e+00    3.286393e-01    7.979069e-02   -2.895422e+00    0.000000e+00    0.000000e+00
            1.000000e+02    6.287160e+01    1.152449e-01    0.000000e+00    0.000000e+00    0.000000e+00   -6.821944e-01    8.195868e-02   -1.628517e+00    0.000000e+00    0.000000e+00
            1.100000e+02    6.287160e+01    1.435680e-01    0.000000e+00    0.000000e+00    0.000000e+00   -1.976259e+00    8.412667e-02   -6.131284e-01    0.000000e+00    0.000000e+00
            1.200000e+02    6.287160e+01    1.718912e-01    0.000000e+00    0.000000e+00    0.000000e+00   -3.553555e+00    8.629465e-02   -2.954865e-01    0.000000e+00    0.000000e+00
        """)


        out = np.loadtxt(output)

        np.testing.assert_array_almost_equal(intF[iE].x[iCase, :], out[:, 0], decimal=3)
        np.testing.assert_array_almost_equal(intF[iE].Nx[iCase, :], out[:, 1], decimal=3)
        np.testing.assert_array_almost_equal(intF[iE].Vy[iCase, :], out[:, 2], decimal=3)
        np.testing.assert_array_almost_equal(intF[iE].Vz[iCase, :], out[:, 3], decimal=3)
        np.testing.assert_array_almost_equal(intF[iE].Tx[iCase, :], out[:, 4], decimal=3)
        np.testing.assert_array_almost_equal(intF[iE].My[iCase, :], out[:, 5], decimal=3)
        np.testing.assert_array_almost_equal(intF[iE].Mz[iCase, :], out[:, 6], decimal=3)
        np.testing.assert_array_almost_equal(intF[iE].Dx[iCase, :], out[:, 7], decimal=3)
        np.testing.assert_array_almost_equal(intF[iE].Dy[iCase, :], out[:, 8], decimal=3)
        np.testing.assert_array_almost_equal(intF[iE].Dz[iCase, :], out[:, 9], decimal=3)
        np.testing.assert_array_almost_equal(intF[iE].Rx[iCase, :], out[:, 10], decimal=3)


    def test_if2(self):

        intF = self.internalForces
        iE = 7
        iCase = 1

        output = StringIO("""
          0.000000e+00   -2.116037e+01   -6.145897e-03    0.000000e+00    0.000000e+00    0.000000e+00    1.708194e-01   -1.067463e+00   -7.196515e-02    0.000000e+00    0.000000e+00
          1.000000e+01   -2.113205e+01   -6.145897e-03    0.000000e+00    0.000000e+00    0.000000e+00    2.322784e-01   -1.068192e+00   -9.033631e-01    0.000000e+00    0.000000e+00
          2.000000e+01   -2.110372e+01   -6.145897e-03    0.000000e+00    0.000000e+00    0.000000e+00    2.937373e-01   -1.068920e+00   -1.654665e+00    0.000000e+00    0.000000e+00
          3.000000e+01   -2.107540e+01   -6.145897e-03    0.000000e+00    0.000000e+00    0.000000e+00    3.551963e-01   -1.069647e+00   -2.304678e+00    0.000000e+00    0.000000e+00
          4.000000e+01   -2.104708e+01   -6.145897e-03    0.000000e+00    0.000000e+00    0.000000e+00    4.166553e-01   -1.070373e+00   -2.832210e+00    0.000000e+00    0.000000e+00
          5.000000e+01   -2.101876e+01   -6.145897e-03    0.000000e+00    0.000000e+00    0.000000e+00    4.781142e-01   -1.071099e+00   -3.216067e+00    0.000000e+00    0.000000e+00
          6.000000e+01   -2.099043e+01   -6.145897e-03    0.000000e+00    0.000000e+00    0.000000e+00    5.395732e-01   -1.071823e+00   -3.435058e+00    0.000000e+00    0.000000e+00
          7.000000e+01   -2.096211e+01   -6.145897e-03    0.000000e+00    0.000000e+00    0.000000e+00    6.010322e-01   -1.072546e+00   -3.467988e+00    0.000000e+00    0.000000e+00
          8.000000e+01   -2.093379e+01   -6.145897e-03    0.000000e+00    0.000000e+00    0.000000e+00    6.624911e-01   -1.073269e+00   -3.293667e+00    0.000000e+00    0.000000e+00
          9.000000e+01   -2.090546e+01   -6.145897e-03    0.000000e+00    0.000000e+00    0.000000e+00    7.239501e-01   -1.073990e+00   -2.890900e+00    0.000000e+00    0.000000e+00
          1.000000e+02   -2.087714e+01   -6.145897e-03    0.000000e+00    0.000000e+00    0.000000e+00    7.854091e-01   -1.074710e+00   -2.238495e+00    0.000000e+00    0.000000e+00
          1.100000e+02   -2.084882e+01   -6.145897e-03    0.000000e+00    0.000000e+00    0.000000e+00    8.468680e-01   -1.075430e+00   -1.315259e+00    0.000000e+00    0.000000e+00
          1.200000e+02   -2.082049e+01   -6.145897e-03    0.000000e+00    0.000000e+00    0.000000e+00    9.083270e-01   -1.076148e+00   -1.000000e-01    0.000000e+00    0.000000e+00
        """)


        out = np.loadtxt(output)

        np.testing.assert_array_almost_equal(intF[iE].x[iCase, :], out[:, 0], decimal=3)
        np.testing.assert_array_almost_equal(intF[iE].Nx[iCase, :], out[:, 1], decimal=3)
        np.testing.assert_array_almost_equal(intF[iE].Vy[iCase, :], out[:, 2], decimal=3)
        np.testing.assert_array_almost_equal(intF[iE].Vz[iCase, :], out[:, 3], decimal=3)
        np.testing.assert_array_almost_equal(intF[iE].Tx[iCase, :], out[:, 4], decimal=3)
        np.testing.assert_array_almost_equal(intF[iE].My[iCase, :], out[:, 5], decimal=3)
        np.testing.assert_array_almost_equal(intF[iE].Mz[iCase, :], out[:, 6], decimal=3)
        np.testing.assert_array_almost_equal(intF[iE].Dx[iCase, :], out[:, 7], decimal=3)
        np.testing.assert_array_almost_equal(intF[iE].Dy[iCase, :], out[:, 8], decimal=3)
        np.testing.assert_array_almost_equal(intF[iE].Dz[iCase, :], out[:, 9], decimal=3)
        np.testing.assert_array_almost_equal(intF[iE].Rx[iCase, :], out[:, 10], decimal=3)






class FrameTestEXB(unittest.TestCase):

    def setUp(self):

        # nodes
        string = StringIO("""
        1   0.0 0.0 1000    0.0
        2   -1200   -900    0.0 0.0
        3    1200   -900    0.0 0.0
        4    1200    900    0.0 0.0
        5   -1200    900    0.0 0.0
        """)
        out = np.loadtxt(string)

        nodes = NodeData(out[:, 0], out[:, 1], out[:, 2], out[:, 3], out[:, 4])

        # reactions
        string = StringIO("""
          2 1  1  1  1  1  1
          3 1  1  1  1  1  1
          4 1  1  1  1  1  1
          5 1  1  1  1  1  1
        """)
        out = np.loadtxt(string, dtype=np.int)

        reactions = ReactionData(out[:, 0], out[:, 1], out[:, 2], out[:, 3], out[:, 4],
            out[:, 5], out[:, 6])

        # elements

        string = StringIO("""
        1 2 1   36.0    20.0    20.0    1000    492     492 200000  79300  0 7.85e-9
        2 1 3   36.0    20.0    20.0    1000    492     492 200000  79300  0 7.85e-9
        3 1 4   36.0    20.0    20.0    1000    492     492 200000  79300  0 7.85e-9
        4 5 1   36.0    20.0    20.0    1000    492     492 200000  79300  0 7.85e-9
        """)
        out = np.loadtxt(string)

        elements = ElementData(out[:, 0], out[:, 1], out[:, 2], out[:, 3], out[:, 4],
            out[:, 5], out[:, 6], out[:, 7], out[:, 8], out[:, 9], out[:, 10],
            out[:, 11], out[:, 12])

        # parameters
        shear = 1               # 1: include shear deformation
        geom = 1                # 1: include geometric stiffness
        exagg_static = 10.0     # exaggerate mesh deformations
        dx = 10.0               # x-axis increment for internal forces
        other = OtherData(shear, geom, exagg_static, dx)


        frame = Frame(nodes, reactions, elements, other)



        # load cases 1
        gx = 0.0
        gy = 0.0
        gz = -9806.33

        load = StaticLoadCase(gx, gy, gz)

        nF = np.array([1])
        Fx = np.array([100.0])
        Fy = np.array([-200.0])
        Fz = np.array([-100.0])
        Mxx = np.array([0.0])
        Myy = np.array([0.0])
        Mzz = np.array([0.0])

        load.changePointLoads(nF, Fx, Fy, Fz, Mxx, Myy, Mzz)

        frame.addLoadCase(load)



        gx = 0.0
        gy = 0.0
        gz = -9806.33

        load = StaticLoadCase(gx, gy, gz)

        EL = np.array([1, 2])
        Px = np.array([0.0, 0.0])
        Py = np.array([100.0, -200.0])
        Pz = np.array([-900.0, 200.0])
        x = np.array([600.0, 800.0])
        load.changeElementLoads(EL, Px, Py, Pz, x)

        frame.addLoadCase(load)





        nM = 6               # number of desired dynamic modes of vibration
        Mmethod = 1                               # 1: subspace Jacobi     2: Stodola
        lump = 0               # 0: consistent mass ... 1: lumped mass matrix
        tol = 1e-9                # mode shape tolerance
        shift = 0.0             # shift value ... for unrestrained structures
        exagg_modal = 10.0                            # exaggerate modal mesh deformations

        dynamic = DynamicAnalysis(nM, Mmethod, lump, tol, shift, exagg_modal)

        N = np.array([1])
        EMs = np.array([0.1])
        EMx = np.array([0.0])
        EMy = np.array([0.0])
        EMz = np.array([0.0])
        EMxy = np.array([0.0])
        EMxz = np.array([0.0])
        EMyz = np.array([0.0])
        rhox = np.array([0.0])
        rhoy = np.array([0.0])
        rhoz = np.array([0.0])
        dynamic.changeExtraInertia(N, EMs, EMx, EMy, EMz, EMxy, EMxz, EMyz, rhox, rhoy, rhoz)

        frame.useDynamicAnalysis(dynamic)



        self.displacements, self.forces, self.reactions, self.internalForces, self.mass, self.modal = frame.run()


    def test_disp1(self):

        disp = self.displacements
        iCase = 0

        node = np.array([1, 2, 3, 4, 5])
        dx = np.array([0.014127, 0.0, 0.0, 0.0, 0.0])
        dy = np.array([-0.050229, 0.0, 0.0, 0.0, 0.0])
        dz = np.array([-0.022374, 0.0, 0.0, 0.0, 0.0])
        dxrot = np.array([0.000037, 0.0, 0.0, 0.0, 0.0])
        dyrot = np.array([0.000009, 0.0, 0.0, 0.0, 0.0])
        dzrot = np.array([0.000000, 0.0, 0.0, 0.0, 0.0])


        np.testing.assert_equal(disp.node[iCase, :], node)
        np.testing.assert_almost_equal(disp.dx[iCase, :], dx, decimal=6)
        np.testing.assert_almost_equal(disp.dy[iCase, :], dy, decimal=6)
        np.testing.assert_almost_equal(disp.dz[iCase, :], dz, decimal=6)
        np.testing.assert_almost_equal(disp.dxrot[iCase, :], dxrot, decimal=6)
        np.testing.assert_almost_equal(disp.dyrot[iCase, :], dyrot, decimal=6)
        np.testing.assert_almost_equal(disp.dzrot[iCase, :], dzrot, decimal=6)



    def test_force1(self):

        forces = self.forces
        iCase = 0

        output = StringIO("""
             1      2    113.543      0.003      2.082     -1.289   -627.688      6.039
             1      1   -110.772     -0.003      2.075      1.289    620.134      4.572
             2      1    185.886     -0.000      2.074      0.904   -620.112     -2.775
             2      3   -188.657      0.000      2.083     -0.904    627.326     -3.504
             3      1    -14.410     -0.007      2.075      1.285   -622.619     -4.567
             3      4     11.639      0.007      2.082     -1.285    628.131     -6.780
             4      5    -86.753      0.006      2.084     -0.908   -629.365      4.619
             4      1     89.524     -0.006      2.073      0.908    623.618      2.765
        """)

        out = np.loadtxt(output)

        np.testing.assert_array_equal(forces.element[iCase, :], out[:, 0])
        np.testing.assert_array_equal(forces.node[iCase, :], out[:, 1])
        np.testing.assert_array_almost_equal(forces.Nx[iCase, :], out[:, 2], decimal=3)
        np.testing.assert_array_almost_equal(forces.Vy[iCase, :], out[:, 3], decimal=3)
        np.testing.assert_array_almost_equal(forces.Vz[iCase, :], out[:, 4], decimal=3)
        np.testing.assert_array_almost_equal(forces.Txx[iCase, :], out[:, 5], decimal=3)
        np.testing.assert_array_almost_equal(forces.Myy[iCase, :], out[:, 6], decimal=3)
        np.testing.assert_array_almost_equal(forces.Mzz[iCase, :], out[:, 7], decimal=3)



    def test_reactions1(self):

        reactions = self.reactions
        iCase = 0

        output = StringIO("""
             1      0.0 0.0 0.0 0.0 0.0 0.0
             2      74.653      55.994      64.715     373.075    -504.804       4.310
             3    -124.653      93.490     106.381     374.239     503.478      -2.414
             4       8.667       6.509      -4.724    -380.743     499.607      -4.929
             5     -58.667      44.008     -46.388    -380.273    -501.502       3.340
        """)

        out = np.loadtxt(output)

        np.testing.assert_array_equal(reactions.node[iCase, :], out[:, 0])
        np.testing.assert_array_almost_equal(reactions.Fx[iCase, :], out[:, 1], decimal=3)
        np.testing.assert_array_almost_equal(reactions.Fy[iCase, :], out[:, 2], decimal=3)
        np.testing.assert_array_almost_equal(reactions.Fz[iCase, :], out[:, 3], decimal=3)
        np.testing.assert_array_almost_equal(reactions.Mxx[iCase, :], out[:, 4], decimal=3)
        np.testing.assert_array_almost_equal(reactions.Myy[iCase, :], out[:, 5], decimal=3)
        np.testing.assert_array_almost_equal(reactions.Mzz[iCase, :], out[:, 6], decimal=3)



    def test_if1(self):

        intF = self.internalForces
        iE = 1
        iCase = 0

        output = StringIO("""
            0.000000e+00     -1.858856e+02    1.892972e-04   -2.073996e+00   -9.040717e-01   -6.201121e+02    2.774643e+00    4.689001e-02   -3.170648e-02    4.370103e-03    2.055282e-05
          1.000000e+01   -1.859010e+02    1.892972e-04   -2.050938e+00   -9.040717e-01   -5.994829e+02    2.739763e+00    4.663183e-02   -3.155635e-02    3.779178e-03    2.043881e-05
          2.000000e+01   -1.859164e+02    1.892972e-04   -2.027879e+00   -9.040717e-01   -5.790842e+02    2.704883e+00    4.637362e-02   -3.140383e-02    2.573358e-03    2.032480e-05
          3.000000e+01   -1.859317e+02    1.892972e-04   -2.004821e+00   -9.040717e-01   -5.589162e+02    2.670003e+00    4.611539e-02   -3.124895e-02    7.799139e-04    2.021080e-05
          4.000000e+01   -1.859471e+02    1.892972e-04   -1.981763e+00   -9.040717e-01   -5.389787e+02    2.635122e+00    4.585714e-02   -3.109174e-02   -1.580660e-03    2.009679e-05
          5.000000e+01   -1.859625e+02    1.892972e-04   -1.958704e+00   -9.040717e-01   -5.192718e+02    2.600242e+00    4.559887e-02   -3.093225e-02   -4.488101e-03    1.998278e-05
          6.000000e+01   -1.859779e+02    1.892972e-04   -1.935646e+00   -9.040717e-01   -4.997955e+02    2.565362e+00    4.534058e-02   -3.077051e-02   -7.922382e-03    1.986878e-05
          7.000000e+01   -1.859932e+02    1.892972e-04   -1.912588e+00   -9.040717e-01   -4.805498e+02    2.530482e+00    4.508227e-02   -3.060655e-02   -1.186371e-02    1.975477e-05
          8.000000e+01   -1.860086e+02    1.892972e-04   -1.889529e+00   -9.040717e-01   -4.615347e+02    2.495601e+00    4.482393e-02   -3.044041e-02   -1.629253e-02    1.964077e-05
          9.000000e+01   -1.860240e+02    1.892972e-04   -1.866471e+00   -9.040717e-01   -4.427501e+02    2.460721e+00    4.456558e-02   -3.027212e-02   -2.118951e-02    1.952676e-05
          1.000000e+02   -1.860394e+02    1.892972e-04   -1.843413e+00   -9.040717e-01   -4.241962e+02    2.425841e+00    4.430720e-02   -3.010173e-02   -2.653556e-02    1.941275e-05
          1.100000e+02   -1.860547e+02    1.892972e-04   -1.820354e+00   -9.040717e-01   -4.058728e+02    2.390961e+00    4.404880e-02   -2.992925e-02   -3.231184e-02    1.929875e-05
          1.200000e+02   -1.860701e+02    1.892972e-04   -1.797296e+00   -9.040717e-01   -3.877800e+02    2.356080e+00    4.379038e-02   -2.975474e-02   -3.849971e-02    1.918474e-05
          1.300000e+02   -1.860855e+02    1.892972e-04   -1.774237e+00   -9.040717e-01   -3.699178e+02    2.321200e+00    4.353194e-02   -2.957823e-02   -4.508079e-02    1.907073e-05
          1.400000e+02   -1.861008e+02    1.892972e-04   -1.751179e+00   -9.040717e-01   -3.522862e+02    2.286320e+00    4.327348e-02   -2.939974e-02   -5.203693e-02    1.895673e-05
          1.500000e+02   -1.861162e+02    1.892972e-04   -1.728121e+00   -9.040717e-01   -3.348851e+02    2.251439e+00    4.301499e-02   -2.921933e-02   -5.935021e-02    1.884272e-05
          1.600000e+02   -1.861316e+02    1.892972e-04   -1.705062e+00   -9.040717e-01   -3.177147e+02    2.216559e+00    4.275649e-02   -2.903701e-02   -6.700294e-02    1.872871e-05
          1.700000e+02   -1.861470e+02    1.892972e-04   -1.682004e+00   -9.040717e-01   -3.007748e+02    2.181679e+00    4.249796e-02   -2.885284e-02   -7.497768e-02    1.861471e-05
          1.800000e+02   -1.861623e+02    1.892972e-04   -1.658946e+00   -9.040717e-01   -2.840655e+02    2.146799e+00    4.223941e-02   -2.866683e-02   -8.325721e-02    1.850070e-05
          1.900000e+02   -1.861777e+02    1.892972e-04   -1.635887e+00   -9.040717e-01   -2.675868e+02    2.111918e+00    4.198084e-02   -2.847904e-02   -9.182455e-02    1.838669e-05
          2.000000e+02   -1.861931e+02    1.892972e-04   -1.612829e+00   -9.040717e-01   -2.513386e+02    2.077038e+00    4.172225e-02   -2.828949e-02   -1.006630e-01    1.827269e-05
          2.100000e+02   -1.862084e+02    1.892972e-04   -1.589771e+00   -9.040717e-01   -2.353211e+02    2.042158e+00    4.146364e-02   -2.809822e-02   -1.097559e-01    1.815868e-05
          2.200000e+02   -1.862238e+02    1.892972e-04   -1.566712e+00   -9.040717e-01   -2.195341e+02    2.007278e+00    4.120501e-02   -2.790526e-02   -1.190871e-01    1.804467e-05
          2.300000e+02   -1.862392e+02    1.892972e-04   -1.543654e+00   -9.040717e-01   -2.039778e+02    1.972397e+00    4.094635e-02   -2.771065e-02   -1.286406e-01    1.793067e-05
          2.400000e+02   -1.862546e+02    1.892972e-04   -1.520596e+00   -9.040717e-01   -1.886520e+02    1.937517e+00    4.068767e-02   -2.751443e-02   -1.384005e-01    1.781666e-05
          2.500000e+02   -1.862699e+02    1.892972e-04   -1.497537e+00   -9.040717e-01   -1.735567e+02    1.902637e+00    4.042898e-02   -2.731663e-02   -1.483512e-01    1.770265e-05
          2.600000e+02   -1.862853e+02    1.892972e-04   -1.474479e+00   -9.040717e-01   -1.586921e+02    1.867757e+00    4.017026e-02   -2.711729e-02   -1.584774e-01    1.758865e-05
          2.700000e+02   -1.863007e+02    1.892972e-04   -1.451421e+00   -9.040717e-01   -1.440581e+02    1.832876e+00    3.991152e-02   -2.691644e-02   -1.687640e-01    1.747464e-05
          2.800000e+02   -1.863161e+02    1.892972e-04   -1.428362e+00   -9.040717e-01   -1.296546e+02    1.797996e+00    3.965276e-02   -2.671412e-02   -1.791962e-01    1.736064e-05
          2.900000e+02   -1.863314e+02    1.892972e-04   -1.405304e+00   -9.040717e-01   -1.154817e+02    1.763116e+00    3.939397e-02   -2.651036e-02   -1.897592e-01    1.724663e-05
          3.000000e+02   -1.863468e+02    1.892972e-04   -1.382246e+00   -9.040717e-01   -1.015394e+02    1.728236e+00    3.913517e-02   -2.630520e-02   -2.004387e-01    1.713262e-05
          3.100000e+02   -1.863622e+02    1.892972e-04   -1.359187e+00   -9.040717e-01   -8.782773e+01    1.693355e+00    3.887634e-02   -2.609868e-02   -2.112205e-01    1.701862e-05
          3.200000e+02   -1.863775e+02    1.892972e-04   -1.336129e+00   -9.040717e-01   -7.434661e+01    1.658475e+00    3.861750e-02   -2.589082e-02   -2.220907e-01    1.690461e-05
          3.300000e+02   -1.863929e+02    1.892972e-04   -1.313070e+00   -9.040717e-01   -6.109606e+01    1.623595e+00    3.835863e-02   -2.568167e-02   -2.330356e-01    1.679060e-05
          3.400000e+02   -1.864083e+02    1.892972e-04   -1.290012e+00   -9.040717e-01   -4.807610e+01    1.588715e+00    3.809974e-02   -2.547126e-02   -2.440417e-01    1.667660e-05
          3.500000e+02   -1.864237e+02    1.892972e-04   -1.266954e+00   -9.040717e-01   -3.528673e+01    1.553834e+00    3.784083e-02   -2.525962e-02   -2.550958e-01    1.656259e-05
          3.600000e+02   -1.864390e+02    1.892972e-04   -1.243895e+00   -9.040717e-01   -2.272794e+01    1.518954e+00    3.758189e-02   -2.504680e-02   -2.661849e-01    1.644858e-05
          3.700000e+02   -1.864544e+02    1.892972e-04   -1.220837e+00   -9.040717e-01   -1.039973e+01    1.484074e+00    3.732294e-02   -2.483282e-02   -2.772962e-01    1.633458e-05
          3.800000e+02   -1.864698e+02    1.892972e-04   -1.197779e+00   -9.040717e-01    1.697898e+00    1.449193e+00    3.706397e-02   -2.461773e-02   -2.884172e-01    1.622057e-05
          3.900000e+02   -1.864852e+02    1.892972e-04   -1.174720e+00   -9.040717e-01    1.356494e+01    1.414313e+00    3.680497e-02   -2.440155e-02   -2.995356e-01    1.610656e-05
          4.000000e+02   -1.865005e+02    1.892972e-04   -1.151662e+00   -9.040717e-01    2.520140e+01    1.379433e+00    3.654595e-02   -2.418433e-02   -3.106393e-01    1.599256e-05
          4.100000e+02   -1.865159e+02    1.892972e-04   -1.128604e+00   -9.040717e-01    3.660727e+01    1.344553e+00    3.628691e-02   -2.396609e-02   -3.217166e-01    1.587855e-05
          4.200000e+02   -1.865313e+02    1.892972e-04   -1.105545e+00   -9.040717e-01    4.778257e+01    1.309672e+00    3.602785e-02   -2.374688e-02   -3.327557e-01    1.576454e-05
          4.300000e+02   -1.865466e+02    1.892972e-04   -1.082487e+00   -9.040717e-01    5.872727e+01    1.274792e+00    3.576877e-02   -2.352673e-02   -3.437455e-01    1.565054e-05
          4.400000e+02   -1.865620e+02    1.892972e-04   -1.059429e+00   -9.040717e-01    6.944140e+01    1.239912e+00    3.550967e-02   -2.330567e-02   -3.546746e-01    1.553653e-05
          4.500000e+02   -1.865774e+02    1.892972e-04   -1.036370e+00   -9.040717e-01    7.992494e+01    1.205032e+00    3.525054e-02   -2.308374e-02   -3.655324e-01    1.542252e-05
          4.600000e+02   -1.865928e+02    1.892972e-04   -1.013312e+00   -9.040717e-01    9.017790e+01    1.170151e+00    3.499140e-02   -2.286098e-02   -3.763080e-01    1.530852e-05
          4.700000e+02   -1.866081e+02    1.892972e-04   -9.902536e-01   -9.040717e-01    1.002003e+02    1.135271e+00    3.473223e-02   -2.263742e-02   -3.869911e-01    1.519451e-05
          4.800000e+02   -1.866235e+02    1.892972e-04   -9.671952e-01   -9.040717e-01    1.099921e+02    1.100391e+00    3.447304e-02   -2.241310e-02   -3.975715e-01    1.508050e-05
          4.900000e+02   -1.866389e+02    1.892972e-04   -9.441369e-01   -9.040717e-01    1.195533e+02    1.065511e+00    3.421383e-02   -2.218805e-02   -4.080392e-01    1.496650e-05
          5.000000e+02   -1.866542e+02    1.892972e-04   -9.210785e-01   -9.040717e-01    1.288839e+02    1.030630e+00    3.395460e-02   -2.196230e-02   -4.183846e-01    1.485249e-05
          5.100000e+02   -1.866696e+02    1.892972e-04   -8.980202e-01   -9.040717e-01    1.379839e+02    9.957501e-01    3.369535e-02   -2.173590e-02   -4.285981e-01    1.473849e-05
          5.200000e+02   -1.866850e+02    1.892972e-04   -8.749618e-01   -9.040717e-01    1.468534e+02    9.608698e-01    3.343607e-02   -2.150888e-02   -4.386706e-01    1.462448e-05
          5.300000e+02   -1.867004e+02    1.892972e-04   -8.519035e-01   -9.040717e-01    1.554923e+02    9.259895e-01    3.317678e-02   -2.128127e-02   -4.485929e-01    1.451047e-05
          5.400000e+02   -1.867157e+02    1.892972e-04   -8.288451e-01   -9.040717e-01    1.639006e+02    8.911093e-01    3.291746e-02   -2.105311e-02   -4.583563e-01    1.439647e-05
          5.500000e+02   -1.867311e+02    1.892972e-04   -8.057868e-01   -9.040717e-01    1.720783e+02    8.562290e-01    3.265812e-02   -2.082444e-02   -4.679522e-01    1.428246e-05
          5.600000e+02   -1.867465e+02    1.892972e-04   -7.827284e-01   -9.040717e-01    1.800254e+02    8.213487e-01    3.239876e-02   -2.059529e-02   -4.773725e-01    1.416845e-05
          5.700000e+02   -1.867619e+02    1.892972e-04   -7.596701e-01   -9.040717e-01    1.877419e+02    7.864685e-01    3.213938e-02   -2.036569e-02   -4.866088e-01    1.405445e-05
          5.800000e+02   -1.867772e+02    1.892972e-04   -7.366117e-01   -9.040717e-01    1.952279e+02    7.515882e-01    3.187998e-02   -2.013568e-02   -4.956536e-01    1.394044e-05
          5.900000e+02   -1.867926e+02    1.892972e-04   -7.135534e-01   -9.040717e-01    2.024832e+02    7.167080e-01    3.162056e-02   -1.990530e-02   -5.044990e-01    1.382643e-05
          6.000000e+02   -1.868080e+02    1.892972e-04   -6.904950e-01   -9.040717e-01    2.095080e+02    6.818277e-01    3.136111e-02   -1.967458e-02   -5.131378e-01    1.371243e-05
          6.100000e+02   -1.868233e+02    1.892972e-04   -6.674367e-01   -9.040717e-01    2.163022e+02    6.469474e-01    3.110164e-02   -1.944356e-02   -5.215628e-01    1.359842e-05
          6.200000e+02   -1.868387e+02    1.892972e-04   -6.443783e-01   -9.040717e-01    2.228659e+02    6.120672e-01    3.084216e-02   -1.921227e-02   -5.297671e-01    1.348441e-05
          6.300000e+02   -1.868541e+02    1.892972e-04   -6.213200e-01   -9.040717e-01    2.291989e+02    5.771869e-01    3.058265e-02   -1.898075e-02   -5.377440e-01    1.337041e-05
          6.400000e+02   -1.868695e+02    1.892972e-04   -5.982616e-01   -9.040717e-01    2.353014e+02    5.423067e-01    3.032312e-02   -1.874904e-02   -5.454872e-01    1.325640e-05
          6.500000e+02   -1.868848e+02    1.892972e-04   -5.752033e-01   -9.040717e-01    2.411732e+02    5.074264e-01    3.006357e-02   -1.851716e-02   -5.529903e-01    1.314239e-05
          6.600000e+02   -1.869002e+02    1.892972e-04   -5.521449e-01   -9.040717e-01    2.468145e+02    4.725461e-01    2.980399e-02   -1.828516e-02   -5.602475e-01    1.302839e-05
          6.700000e+02   -1.869156e+02    1.892972e-04   -5.290866e-01   -9.040717e-01    2.522252e+02    4.376659e-01    2.954440e-02   -1.805306e-02   -5.672529e-01    1.291438e-05
          6.800000e+02   -1.869309e+02    1.892972e-04   -5.060282e-01   -9.040717e-01    2.574053e+02    4.027856e-01    2.928478e-02   -1.782092e-02   -5.740012e-01    1.280037e-05
          6.900000e+02   -1.869463e+02    1.892972e-04   -4.829699e-01   -9.040717e-01    2.623549e+02    3.679053e-01    2.902515e-02   -1.758875e-02   -5.804870e-01    1.268637e-05
          7.000000e+02   -1.869617e+02    1.892972e-04   -4.599115e-01   -9.040717e-01    2.670738e+02    3.330251e-01    2.876549e-02   -1.735660e-02   -5.867053e-01    1.257236e-05
          7.100000e+02   -1.869771e+02    1.892972e-04   -4.368532e-01   -9.040717e-01    2.715622e+02    2.981448e-01    2.850581e-02   -1.712451e-02   -5.926513e-01    1.245835e-05
          7.200000e+02   -1.869924e+02    1.892972e-04   -4.137948e-01   -9.040717e-01    2.758200e+02    2.632646e-01    2.824611e-02   -1.689250e-02   -5.983205e-01    1.234435e-05
          7.300000e+02   -1.870078e+02    1.892972e-04   -3.907365e-01   -9.040717e-01    2.798472e+02    2.283843e-01    2.798639e-02   -1.666061e-02   -6.037084e-01    1.223034e-05
          7.400000e+02   -1.870232e+02    1.892972e-04   -3.676781e-01   -9.040717e-01    2.836438e+02    1.935040e-01    2.772664e-02   -1.642888e-02   -6.088111e-01    1.211634e-05
          7.500000e+02   -1.870386e+02    1.892972e-04   -3.446198e-01   -9.040717e-01    2.872098e+02    1.586238e-01    2.746688e-02   -1.619735e-02   -6.136247e-01    1.200233e-05
          7.600000e+02   -1.870539e+02    1.892972e-04   -3.215614e-01   -9.040717e-01    2.905453e+02    1.237435e-01    2.720709e-02   -1.596604e-02   -6.181455e-01    1.188832e-05
          7.700000e+02   -1.870693e+02    1.892972e-04   -2.985031e-01   -9.040717e-01    2.936502e+02    8.886324e-02    2.694728e-02   -1.573500e-02   -6.223702e-01    1.177432e-05
          7.800000e+02   -1.870847e+02    1.892972e-04   -2.754447e-01   -9.040717e-01    2.965245e+02    5.398298e-02    2.668745e-02   -1.550426e-02   -6.262956e-01    1.166031e-05
          7.900000e+02   -1.871000e+02    1.892972e-04   -2.523864e-01   -9.040717e-01    2.991682e+02    1.910272e-02    2.642760e-02   -1.527386e-02   -6.299187e-01    1.154630e-05
          8.000000e+02   -1.871154e+02    1.892972e-04   -2.293280e-01   -9.040717e-01    3.015813e+02   -1.577754e-02    2.616773e-02   -1.504382e-02   -6.332370e-01    1.143230e-05
          8.100000e+02   -1.871308e+02    1.892972e-04   -2.062697e-01   -9.040717e-01    3.037638e+02   -5.065780e-02    2.590784e-02   -1.481420e-02   -6.362479e-01    1.131829e-05
          8.200000e+02   -1.871462e+02    1.892972e-04   -1.832113e-01   -9.040717e-01    3.057158e+02   -8.553807e-02    2.564792e-02   -1.458501e-02   -6.389492e-01    1.120428e-05
          8.300000e+02   -1.871615e+02    1.892972e-04   -1.601530e-01   -9.040717e-01    3.074371e+02   -1.204183e-01    2.538799e-02   -1.435630e-02   -6.413389e-01    1.109028e-05
          8.400000e+02   -1.871769e+02    1.892972e-04   -1.370946e-01   -9.040717e-01    3.089279e+02   -1.552986e-01    2.512803e-02   -1.412811e-02   -6.434153e-01    1.097627e-05
          8.500000e+02   -1.871923e+02    1.892972e-04   -1.140363e-01   -9.040717e-01    3.101881e+02   -1.901789e-01    2.486805e-02   -1.390046e-02   -6.451770e-01    1.086226e-05
          8.600000e+02   -1.872076e+02    1.892972e-04   -9.097793e-02   -9.040717e-01    3.112177e+02   -2.250591e-01    2.460805e-02   -1.367340e-02   -6.466225e-01    1.074826e-05
          8.700000e+02   -1.872230e+02    1.892972e-04   -6.791958e-02   -9.040717e-01    3.120168e+02   -2.599394e-01    2.434803e-02   -1.344695e-02   -6.477508e-01    1.063425e-05
          8.800000e+02   -1.872384e+02    1.892972e-04   -4.486123e-02   -9.040717e-01    3.125852e+02   -2.948196e-01    2.408799e-02   -1.322116e-02   -6.485612e-01    1.052024e-05
          8.900000e+02   -1.872538e+02    1.892972e-04   -2.180288e-02   -9.040717e-01    3.129231e+02   -3.296999e-01    2.382792e-02   -1.299606e-02   -6.490530e-01    1.040624e-05
          9.000000e+02   -1.872691e+02    1.892972e-04    1.255467e-03   -9.040717e-01    3.130304e+02   -3.645802e-01    2.356784e-02   -1.277169e-02   -6.492260e-01    1.029223e-05
          9.100000e+02   -1.872845e+02    1.892972e-04    2.431382e-02   -9.040717e-01    3.129071e+02   -3.994604e-01    2.330773e-02   -1.254807e-02   -6.490799e-01    1.017822e-05
          9.200000e+02   -1.872999e+02    1.892972e-04    4.737217e-02   -9.040717e-01    3.125532e+02   -4.343407e-01    2.304760e-02   -1.232526e-02   -6.486150e-01    1.006422e-05
          9.300000e+02   -1.873153e+02    1.892972e-04    7.043052e-02   -9.040717e-01    3.119687e+02   -4.692209e-01    2.278745e-02   -1.210327e-02   -6.478316e-01    9.950211e-06
          9.400000e+02   -1.873306e+02    1.892972e-04    9.348886e-02   -9.040717e-01    3.111537e+02   -5.041012e-01    2.252728e-02   -1.188215e-02   -6.467303e-01    9.836205e-06
          9.500000e+02   -1.873460e+02    1.892972e-04    1.165472e-01   -9.040717e-01    3.101080e+02   -5.389815e-01    2.226709e-02   -1.166194e-02   -6.453118e-01    9.722198e-06
          9.600000e+02   -1.873614e+02    1.892972e-04    1.396056e-01   -9.040717e-01    3.088318e+02   -5.738617e-01    2.200688e-02   -1.144266e-02   -6.435774e-01    9.608192e-06
          9.700000e+02   -1.873767e+02    1.892972e-04    1.626639e-01   -9.040717e-01    3.073250e+02   -6.087420e-01    2.174664e-02   -1.122435e-02   -6.415282e-01    9.494185e-06
          9.800000e+02   -1.873921e+02    1.892972e-04    1.857223e-01   -9.040717e-01    3.055876e+02   -6.436223e-01    2.148639e-02   -1.100706e-02   -6.391658e-01    9.380179e-06
          9.900000e+02   -1.874075e+02    1.892972e-04    2.087806e-01   -9.040717e-01    3.036197e+02   -6.785025e-01    2.122611e-02   -1.079081e-02   -6.364920e-01    9.266172e-06
          1.000000e+03   -1.874229e+02    1.892972e-04    2.318390e-01   -9.040717e-01    3.014211e+02   -7.133828e-01    2.096581e-02   -1.057564e-02   -6.335088e-01    9.152166e-06
          1.010000e+03   -1.874382e+02    1.892972e-04    2.548973e-01   -9.040717e-01    2.989920e+02   -7.482630e-01    2.070549e-02   -1.036158e-02   -6.302183e-01    9.038159e-06
          1.020000e+03   -1.874536e+02    1.892972e-04    2.779557e-01   -9.040717e-01    2.963323e+02   -7.831433e-01    2.044515e-02   -1.014868e-02   -6.266231e-01    8.924153e-06
          1.030000e+03   -1.874690e+02    1.892972e-04    3.010140e-01   -9.040717e-01    2.934420e+02   -8.180236e-01    2.018479e-02   -9.936962e-03   -6.227259e-01    8.810146e-06
          1.040000e+03   -1.874843e+02    1.892972e-04    3.240724e-01   -9.040717e-01    2.903211e+02   -8.529038e-01    1.992440e-02   -9.726465e-03   -6.185297e-01    8.696140e-06
          1.050000e+03   -1.874997e+02    1.892972e-04    3.471307e-01   -9.040717e-01    2.869696e+02   -8.877841e-01    1.966400e-02   -9.517226e-03   -6.140375e-01    8.582133e-06
          1.060000e+03   -1.875151e+02    1.892972e-04    3.701891e-01   -9.040717e-01    2.833876e+02   -9.226644e-01    1.940357e-02   -9.309280e-03   -6.092528e-01    8.468127e-06
          1.070000e+03   -1.875305e+02    1.892972e-04    3.932474e-01   -9.040717e-01    2.795749e+02   -9.575446e-01    1.914312e-02   -9.102661e-03   -6.041792e-01    8.354120e-06
          1.080000e+03   -1.875458e+02    1.892972e-04    4.163058e-01   -9.040717e-01    2.755317e+02   -9.924249e-01    1.888265e-02   -8.897407e-03   -5.988206e-01    8.240113e-06
          1.090000e+03   -1.875612e+02    1.892972e-04    4.393641e-01   -9.040717e-01    2.712579e+02   -1.027305e+00    1.862216e-02   -8.693551e-03   -5.931811e-01    8.126107e-06
          1.100000e+03   -1.875766e+02    1.892972e-04    4.624225e-01   -9.040717e-01    2.667535e+02   -1.062185e+00    1.836165e-02   -8.491129e-03   -5.872651e-01    8.012100e-06
          1.110000e+03   -1.875920e+02    1.892972e-04    4.854808e-01   -9.040717e-01    2.620185e+02   -1.097066e+00    1.810111e-02   -8.290177e-03   -5.810771e-01    7.898094e-06
          1.120000e+03   -1.876073e+02    1.892972e-04    5.085392e-01   -9.040717e-01    2.570530e+02   -1.131946e+00    1.784056e-02   -8.090731e-03   -5.746220e-01    7.784087e-06
          1.130000e+03   -1.876227e+02    1.892972e-04    5.315975e-01   -9.040717e-01    2.518568e+02   -1.166826e+00    1.757998e-02   -7.892825e-03   -5.679048e-01    7.670081e-06
          1.140000e+03   -1.876381e+02    1.892972e-04    5.546559e-01   -9.040717e-01    2.464301e+02   -1.201706e+00    1.731938e-02   -7.696496e-03   -5.609307e-01    7.556074e-06
          1.150000e+03   -1.876534e+02    1.892972e-04    5.777142e-01   -9.040717e-01    2.407728e+02   -1.236587e+00    1.705877e-02   -7.501778e-03   -5.537053e-01    7.442068e-06
          1.160000e+03   -1.876688e+02    1.892972e-04    6.007726e-01   -9.040717e-01    2.348849e+02   -1.271467e+00    1.679813e-02   -7.308707e-03   -5.462344e-01    7.328061e-06
          1.170000e+03   -1.876842e+02    1.892972e-04    6.238309e-01   -9.040717e-01    2.287665e+02   -1.306347e+00    1.653746e-02   -7.117319e-03   -5.385239e-01    7.214055e-06
          1.180000e+03   -1.876996e+02    1.892972e-04    6.468893e-01   -9.040717e-01    2.224174e+02   -1.341228e+00    1.627678e-02   -6.927648e-03   -5.305800e-01    7.100048e-06
          1.190000e+03   -1.877149e+02    1.892972e-04    6.699476e-01   -9.040717e-01    2.158378e+02   -1.376108e+00    1.601608e-02   -6.739732e-03   -5.224092e-01    6.986042e-06
          1.200000e+03   -1.877303e+02    1.892972e-04    6.930060e-01   -9.040717e-01    2.090276e+02   -1.410988e+00    1.575535e-02   -6.553604e-03   -5.140182e-01    6.872035e-06
          1.210000e+03   -1.877457e+02    1.892972e-04    7.160643e-01   -9.040717e-01    2.019867e+02   -1.445868e+00    1.549460e-02   -6.369300e-03   -5.054139e-01    6.758029e-06
          1.220000e+03   -1.877610e+02    1.892972e-04    7.391227e-01   -9.040717e-01    1.947154e+02   -1.480749e+00    1.523383e-02   -6.186856e-03   -4.966034e-01    6.644022e-06
          1.230000e+03   -1.877764e+02    1.892972e-04    7.621810e-01   -9.040717e-01    1.872134e+02   -1.515629e+00    1.497304e-02   -6.006308e-03   -4.875942e-01    6.530016e-06
          1.240000e+03   -1.877918e+02    1.892972e-04    7.852394e-01   -9.040717e-01    1.794808e+02   -1.550509e+00    1.471223e-02   -5.827690e-03   -4.783939e-01    6.416009e-06
          1.250000e+03   -1.878072e+02    1.892972e-04    8.082977e-01   -9.040717e-01    1.715177e+02   -1.585389e+00    1.445140e-02   -5.651038e-03   -4.690102e-01    6.302003e-06
          1.260000e+03   -1.878225e+02    1.892972e-04    8.313561e-01   -9.040717e-01    1.633240e+02   -1.620270e+00    1.419055e-02   -5.476387e-03   -4.594514e-01    6.187996e-06
          1.270000e+03   -1.878379e+02    1.892972e-04    8.544144e-01   -9.040717e-01    1.548997e+02   -1.655150e+00    1.392967e-02   -5.303774e-03   -4.497258e-01    6.073990e-06
          1.280000e+03   -1.878533e+02    1.892972e-04    8.774728e-01   -9.040717e-01    1.462448e+02   -1.690030e+00    1.366877e-02   -5.133233e-03   -4.398418e-01    5.959983e-06
          1.290000e+03   -1.878687e+02    1.892972e-04    9.005311e-01   -9.040717e-01    1.373593e+02   -1.724910e+00    1.340786e-02   -4.964800e-03   -4.298084e-01    5.845977e-06
          1.300000e+03   -1.878840e+02    1.892972e-04    9.235895e-01   -9.040717e-01    1.282432e+02   -1.759791e+00    1.314692e-02   -4.798511e-03   -4.196345e-01    5.731970e-06
          1.310000e+03   -1.878994e+02    1.892972e-04    9.466478e-01   -9.040717e-01    1.188966e+02   -1.794671e+00    1.288596e-02   -4.634400e-03   -4.093293e-01    5.617964e-06
          1.320000e+03   -1.879148e+02    1.892972e-04    9.697061e-01   -9.040717e-01    1.093194e+02   -1.829551e+00    1.262497e-02   -4.472503e-03   -3.989025e-01    5.503957e-06
          1.330000e+03   -1.879301e+02    1.892972e-04    9.927645e-01   -9.040717e-01    9.951158e+01   -1.864431e+00    1.236397e-02   -4.312857e-03   -3.883637e-01    5.389950e-06
          1.340000e+03   -1.879455e+02    1.892972e-04    1.015823e+00   -9.040717e-01    8.947319e+01   -1.899312e+00    1.210295e-02   -4.155495e-03   -3.777229e-01    5.275944e-06
          1.350000e+03   -1.879609e+02    1.892972e-04    1.038881e+00   -9.040717e-01    7.920422e+01   -1.934192e+00    1.184190e-02   -4.000454e-03   -3.669903e-01    5.161937e-06
          1.360000e+03   -1.879763e+02    1.892972e-04    1.061940e+00   -9.040717e-01    6.870466e+01   -1.969072e+00    1.158083e-02   -3.847769e-03   -3.561763e-01    5.047931e-06
          1.370000e+03   -1.879916e+02    1.892972e-04    1.084998e+00   -9.040717e-01    5.797452e+01   -2.003952e+00    1.131974e-02   -3.697475e-03   -3.452917e-01    4.933924e-06
          1.380000e+03   -1.880070e+02    1.892972e-04    1.108056e+00   -9.040717e-01    4.701379e+01   -2.038833e+00    1.105863e-02   -3.549608e-03   -3.343472e-01    4.819918e-06
          1.390000e+03   -1.880224e+02    1.892972e-04    1.131115e+00   -9.040717e-01    3.582249e+01   -2.073713e+00    1.079750e-02   -3.404204e-03   -3.233541e-01    4.705911e-06
          1.400000e+03   -1.880377e+02    1.892972e-04    1.154173e+00   -9.040717e-01    2.440060e+01   -2.108593e+00    1.053635e-02   -3.261298e-03   -3.123237e-01    4.591905e-06
          1.410000e+03   -1.880531e+02    1.892972e-04    1.177231e+00   -9.040717e-01    1.274812e+01   -2.143474e+00    1.027517e-02   -3.120924e-03   -3.012676e-01    4.477898e-06
          1.420000e+03   -1.880685e+02    1.892972e-04    1.200290e+00   -9.040717e-01    8.650633e-01   -2.178354e+00    1.001398e-02   -2.983120e-03   -2.901977e-01    4.363892e-06
          1.430000e+03   -1.880839e+02    1.892972e-04    1.223348e+00   -9.040717e-01   -1.124858e+01   -2.213234e+00    9.752762e-03   -2.847920e-03   -2.791261e-01    4.249885e-06
          1.440000e+03   -1.880992e+02    1.892972e-04    1.246406e+00   -9.040717e-01   -2.359280e+01   -2.248114e+00    9.491524e-03   -2.715359e-03   -2.680650e-01    4.135879e-06
          1.450000e+03   -1.881146e+02    1.892972e-04    1.269465e+00   -9.040717e-01   -3.616761e+01   -2.282995e+00    9.230265e-03   -2.585473e-03   -2.570270e-01    4.021872e-06
          1.460000e+03   -1.881300e+02    1.892972e-04    1.292523e+00   -9.040717e-01   -4.897300e+01   -2.317875e+00    8.968984e-03   -2.458298e-03   -2.460249e-01    3.907866e-06
          1.470000e+03   -1.881454e+02    1.892972e-04    1.315581e+00   -9.040717e-01   -6.200898e+01   -2.352755e+00    8.707681e-03   -2.333869e-03   -2.350716e-01    3.793859e-06
          1.480000e+03   -1.881607e+02    1.892972e-04    1.338640e+00   -9.040717e-01   -7.527554e+01   -2.387635e+00    8.446358e-03   -2.212221e-03   -2.241806e-01    3.679853e-06
          1.490000e+03   -1.881761e+02    1.892972e-04    1.361698e+00   -9.040717e-01   -8.877268e+01   -2.422516e+00    8.185013e-03   -2.093390e-03   -2.133651e-01    3.565846e-06
          1.500000e+03   -1.881915e+02    1.892972e-04    1.384756e+00   -9.040717e-01   -1.025004e+02   -2.457396e+00    7.923646e-03   -1.977411e-03   -2.026390e-01    3.451840e-06
          1.510000e+03   -1.882068e+02    1.892972e-04    1.407815e+00   -9.040717e-01   -1.164587e+02   -2.492276e+00    7.662259e-03   -1.864321e-03   -1.920162e-01    3.337833e-06
          1.520000e+03   -1.882222e+02    1.892972e-04    1.430873e+00   -9.040717e-01   -1.306476e+02   -2.527156e+00    7.400849e-03   -1.754153e-03   -1.815108e-01    3.223827e-06
          1.530000e+03   -1.882376e+02    1.892972e-04    1.453931e+00   -9.040717e-01   -1.450671e+02   -2.562037e+00    7.139419e-03   -1.646944e-03   -1.711374e-01    3.109820e-06
          1.540000e+03   -1.882530e+02    1.892972e-04    1.476990e+00   -9.040717e-01   -1.597171e+02   -2.596917e+00    6.877967e-03   -1.542729e-03   -1.609105e-01    2.995814e-06
          1.550000e+03   -1.882683e+02    1.892972e-04    1.500048e+00   -9.040717e-01   -1.745978e+02   -2.631797e+00    6.616494e-03   -1.441543e-03   -1.508450e-01    2.881807e-06
          1.560000e+03   -1.882837e+02    1.892972e-04    1.523107e+00   -9.040717e-01   -1.897090e+02   -2.666677e+00    6.355000e-03   -1.343423e-03   -1.409561e-01    2.767801e-06
          1.570000e+03   -1.882991e+02    1.892972e-04    1.546165e+00   -9.040717e-01   -2.050508e+02   -2.701558e+00    6.093484e-03   -1.248403e-03   -1.312592e-01    2.653794e-06
          1.580000e+03   -1.883144e+02    1.892972e-04    1.569223e+00   -9.040717e-01   -2.206232e+02   -2.736438e+00    5.831947e-03   -1.156519e-03   -1.217697e-01    2.539788e-06
          1.590000e+03   -1.883298e+02    1.892972e-04    1.592282e+00   -9.040717e-01   -2.364262e+02   -2.771318e+00    5.570388e-03   -1.067806e-03   -1.125036e-01    2.425781e-06
          1.600000e+03   -1.883452e+02    1.892972e-04    1.615340e+00   -9.040717e-01   -2.524598e+02   -2.806199e+00    5.308808e-03   -9.823000e-04   -1.034768e-01    2.311774e-06
          1.610000e+03   -1.883606e+02    1.892972e-04    1.638398e+00   -9.040717e-01   -2.687239e+02   -2.841079e+00    5.047207e-03   -9.000362e-04   -9.470577e-02    2.197768e-06
          1.620000e+03   -1.883759e+02    1.892972e-04    1.661457e+00   -9.040717e-01   -2.852186e+02   -2.875959e+00    4.785585e-03   -8.210501e-04   -8.620695e-02    2.083761e-06
          1.630000e+03   -1.883913e+02    1.892972e-04    1.684515e+00   -9.040717e-01   -3.019439e+02   -2.910839e+00    4.523941e-03   -7.453770e-04   -7.799710e-02    1.969755e-06
          1.640000e+03   -1.884067e+02    1.892972e-04    1.707573e+00   -9.040717e-01   -3.188998e+02   -2.945720e+00    4.262275e-03   -6.730525e-04   -7.009323e-02    1.855748e-06
          1.650000e+03   -1.884221e+02    1.892972e-04    1.730632e+00   -9.040717e-01   -3.360863e+02   -2.980600e+00    4.000589e-03   -6.041121e-04   -6.251257e-02    1.741742e-06
          1.660000e+03   -1.884374e+02    1.892972e-04    1.753690e+00   -9.040717e-01   -3.535034e+02   -3.015480e+00    3.738881e-03   -5.385911e-04   -5.527259e-02    1.627735e-06
          1.670000e+03   -1.884528e+02    1.892972e-04    1.776748e+00   -9.040717e-01   -3.711510e+02   -3.050360e+00    3.477151e-03   -4.765250e-04   -4.839098e-02    1.513729e-06
          1.680000e+03   -1.884682e+02    1.892972e-04    1.799807e+00   -9.040717e-01   -3.890293e+02   -3.085241e+00    3.215401e-03   -4.179492e-04   -4.188569e-02    1.399722e-06
          1.690000e+03   -1.884835e+02    1.892972e-04    1.822865e+00   -9.040717e-01   -4.071381e+02   -3.120121e+00    2.953629e-03   -3.628993e-04   -3.577487e-02    1.285716e-06
          1.700000e+03   -1.884989e+02    1.892972e-04    1.845923e+00   -9.040717e-01   -4.254775e+02   -3.155001e+00    2.691835e-03   -3.114106e-04   -3.007694e-02    1.171709e-06
          1.710000e+03   -1.885143e+02    1.892972e-04    1.868982e+00   -9.040717e-01   -4.440474e+02   -3.189881e+00    2.430021e-03   -2.635186e-04   -2.481053e-02    1.057703e-06
          1.720000e+03   -1.885297e+02    1.892972e-04    1.892040e+00   -9.040717e-01   -4.628480e+02   -3.224762e+00    2.168185e-03   -2.192587e-04   -1.999451e-02    9.436962e-07
          1.730000e+03   -1.885450e+02    1.892972e-04    1.915098e+00   -9.040717e-01   -4.818792e+02   -3.259642e+00    1.906327e-03   -1.786664e-04   -1.564799e-02    8.296897e-07
          1.740000e+03   -1.885604e+02    1.892972e-04    1.938157e+00   -9.040717e-01   -5.011409e+02   -3.294522e+00    1.644448e-03   -1.417772e-04   -1.179031e-02    7.156832e-07
          1.750000e+03   -1.885758e+02    1.892972e-04    1.961215e+00   -9.040717e-01   -5.206332e+02   -3.329402e+00    1.382548e-03   -1.086264e-04   -8.441047e-03    6.016767e-07
          1.760000e+03   -1.885911e+02    1.892972e-04    1.984274e+00   -9.040717e-01   -5.403561e+02   -3.364283e+00    1.120627e-03   -7.924958e-05   -5.620005e-03    4.876702e-07
          1.770000e+03   -1.886065e+02    1.892972e-04    2.007332e+00   -9.040717e-01   -5.603096e+02   -3.399163e+00    8.586840e-04   -5.368214e-05   -3.347230e-03    3.736636e-07
          1.780000e+03   -1.886219e+02    1.892972e-04    2.030390e+00   -9.040717e-01   -5.804936e+02   -3.434043e+00    5.967198e-04   -3.195951e-05   -1.643000e-03    2.596571e-07
          1.790000e+03   -1.886373e+02    1.892972e-04    2.053449e+00   -9.040717e-01   -6.009083e+02   -3.468923e+00    3.347343e-04   -1.411717e-05   -5.278276e-04    1.456506e-07
          1.802776e+03   -1.886569e+02    1.892972e-04    2.082907e+00   -9.040717e-01   -6.273260e+02   -3.504329e+00    0.000000e+00    0.000000e+00    4.336809e-19    0.000000e+00
        """)


        out = np.loadtxt(output)

        np.testing.assert_array_almost_equal(intF[iE].x[iCase, :], out[:, 0], decimal=3)
        np.testing.assert_array_almost_equal(intF[iE].Nx[iCase, :], out[:, 1], decimal=3)
        np.testing.assert_array_almost_equal(intF[iE].Vy[iCase, :], out[:, 2], decimal=3)
        np.testing.assert_array_almost_equal(intF[iE].Vz[iCase, :], out[:, 3], decimal=3)
        np.testing.assert_array_almost_equal(intF[iE].Tx[iCase, :], out[:, 4], decimal=3)
        np.testing.assert_array_almost_equal(intF[iE].My[iCase, :], out[:, 5], decimal=3)
        np.testing.assert_array_almost_equal(intF[iE].Mz[iCase, :], out[:, 6], decimal=3)
        np.testing.assert_array_almost_equal(intF[iE].Dx[iCase, :], out[:, 7], decimal=3)
        np.testing.assert_array_almost_equal(intF[iE].Dy[iCase, :], out[:, 8], decimal=3)
        np.testing.assert_array_almost_equal(intF[iE].Dz[iCase, :], out[:, 9], decimal=3)
        np.testing.assert_array_almost_equal(intF[iE].Rx[iCase, :], out[:, 10], decimal=3)



    def test_mass(self):

        mass = self.mass

        string = StringIO("""
         1 1.00723e-01 1.00738e-01 1.00733e-01 3.51392e+01 4.73634e+01 4.36768e+01
         2 1.26482e+02 1.26482e+02 1.26482e+02 1.26482e+02 1.26482e+02 1.26482e+02
         3 1.26482e+02 1.26482e+02 1.26482e+02 1.26482e+02 1.26482e+02 1.26482e+02
         4 1.26482e+02 1.26482e+02 1.26482e+02 1.26482e+02 1.26482e+02 1.26482e+02
         5 1.26482e+02 1.26482e+02 1.26482e+02 1.26482e+02 1.26482e+02 1.26482e+02
        """)
        out = np.loadtxt(string)

        np.testing.assert_almost_equal(mass.total_mass, 1.020379e-01, decimal=6)
        np.testing.assert_almost_equal(mass.struct_mass, 2.037858e-03, decimal=6)
        np.testing.assert_array_equal(mass.node, out[:, 0])
        np.testing.assert_array_almost_equal(mass.xmass, out[:, 1], decimal=3)
        np.testing.assert_array_almost_equal(mass.ymass, out[:, 2], decimal=3)
        np.testing.assert_array_almost_equal(mass.zmass, out[:, 3], decimal=3)
        np.testing.assert_array_almost_equal(mass.xinrta, out[:, 4], decimal=3)
        np.testing.assert_array_almost_equal(mass.yinrta, out[:, 5], decimal=3)
        np.testing.assert_array_almost_equal(mass.zinrta, out[:, 6], decimal=3)



    def test_modal(self):

        modal = self.modal
        iM = 0

        string = StringIO("""
         1  -2.609e-02  -1.255e-04  -6.331e-06  -2.851e-04   1.344e-01  -5.727e-02
         2  -5.705e-09  -7.082e-10   6.815e-10   4.135e-07  -1.724e-06   1.206e-06
         3  -5.685e-09   6.844e-10  -6.879e-10  -4.094e-07  -1.721e-06   1.210e-06
         4  -3.947e-09  -3.078e-09  -7.076e-10   1.142e-06  -1.155e-06  -9.324e-08
         5  -3.931e-09   3.055e-09   7.135e-10  -1.138e-06  -1.151e-06  -9.690e-08
        """)
        out = np.loadtxt(string)


        np.testing.assert_almost_equal(modal.freq[iM], 18.808236, decimal=5)
        np.testing.assert_almost_equal(modal.xmpf[iM], -2.5467e-02, decimal=6)
        np.testing.assert_almost_equal(modal.ympf[iM], -6.1093e-05, decimal=9)
        np.testing.assert_almost_equal(modal.zmpf[iM], -6.3967e-07, decimal=11)
        np.testing.assert_array_equal(modal.node[iM, :], out[:, 0])
        np.testing.assert_array_almost_equal(modal.xdsp[iM, :], out[:, 1], decimal=3)
        np.testing.assert_array_almost_equal(modal.ydsp[iM, :], out[:, 2], decimal=3)
        np.testing.assert_array_almost_equal(modal.zdsp[iM, :], out[:, 3], decimal=3)
        np.testing.assert_array_almost_equal(modal.xrot[iM, :], out[:, 4], decimal=3)
        np.testing.assert_array_almost_equal(modal.yrot[iM, :], out[:, 5], decimal=3)
        np.testing.assert_array_almost_equal(modal.zrot[iM, :], out[:, 6], decimal=3)





if __name__ == "__main__":
    unittest.main()

    # from unittest import TestSuite, TextTestRunner
    # f = TestSuite()
    # f.addTest(FrameTestEXB('test_modal'))
    # TextTestRunner().run(f)




