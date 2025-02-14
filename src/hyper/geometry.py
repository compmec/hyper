import numpy as np

########################################################################
#                           Gauss Quadrature                           #
########################################################################
#
#
# In 1D, we have: ----- Integration in the interval [-1, 1]
#    ________________________________________________________________
#   | Nb Points |    Position of the points   |         Weight       |
#   |-----------|-----------------------------|----------------------|
#   |     1     |              0              |           2          |
#   |-----------|-----------------------------|----------------------|
#   |     2     |           -a,  a            |         1,  1        |
#   |-----------|-----------------------------|----------------------|
#   |     3     |          -b, 0, b           |    5/9, 8/9, 5/9     |
#    ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾
# Where a = 1/sqrt(3) and b = sqrt(3/5)
#
#
# In 2D, we have: ----- Integration in the interval [-1, 1] x [-1, 1]
#    ________________________________________________________________
#   | Nb Points |    Position of the points   |         Weight       |
#   |-----------|-----------------------------|----------------------|
#   |     1     |           (0, 0)            |           4          |
#   |-----------|-----------------------------|----------------------|
#   |     4     |       (-a, a), (a, a),      |         1, 1,        |
#   |           |      (-a, -a), (a, -a)      |         1, 1         |
#   |-----------|-----------------------------|----------------------|
#   |           |   (-b, b), (0, b), (b, b)   | 25/81, 40/81, 25/81  |
#   |     9     |   (-b, 0), (0, 0), (b, 0)   | 40/81, 64/81, 40/81  |
#   |           |  (-b, -b), (0, -b), (b, -b) | 25/81, 40/81, 25/81  |
#    ‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾‾


class IPTri:
    """
    Numerical integration rules for triangles
    Vertices:
        P1 = (0, 0)
        P2 = (1, 0)
        P3 = (0, 1)
    """

    @staticmethod
    def getX(npg=3):
        if npg == 1:
            X = [(1.0 / 3.0, 1 / 3.0)]
        elif npg == 3:
            X = [(2 / 3, 1 / 6), (1 / 6, 1 / 6), (1 / 6, 2 / 3)]
        elif npg == 4:
            X = [
                (1 / 3, 1 / 3),
                (3 / 5, 1 / 5),
                (1 / 5, 1 / 5),
                (1 / 5, 3 / 5),
            ]
        else:
            raise ValueError("Not found npg = " + str(npg))
        return X

    @staticmethod
    def getW(npg=3):
        if npg == 1:
            W = [1.0 / 2.0]
        elif npg == 3:
            W = [1 / 6, 1 / 6, 1 / 6]
        elif npg == 4:
            W = [-9 / 32, 25 / 96, 25 / 96, 25 / 96]
        else:
            raise ValueError("Not found npg = " + str(npg))
        return W


class IPQua:
    """
    Numerical integration rules for quadrangles
    Vertices:
        P1 = (-1, -1)
        P2 = ( 1, -1)
        P3 = ( 1,  1)
        P4 = (-1,  1)
    """

    @staticmethod
    def getX(npg=3):
        if npg == 1:
            X = [(0, 0)]
        elif npg == 4:
            a = 1 / np.sqrt(3)
            X = [(-a, a), (-a, -a), (a, -a), (a, a)]
        elif npg == 9:
            b = np.sqrt(3 / 5)
            X = [
                (-b, b),
                (0, b),
                (b, b),
                (-b, 0),
                (0, 0),
                (b, 0),
                (-b, -b),
                (0, -b),
                (b, -b),
            ]
        else:
            raise ValueError("Not found npg = " + str(npg))
        return X

    @staticmethod
    def getW(npg=4):
        if npg == 1:
            W = [4]
        elif npg == 4:
            W = [1, 1, 1, 1]
        elif npg == 9:
            W = [
                25 / 81,
                40 / 81,
                25 / 81,
                40 / 81,
                64 / 81,
                40 / 81,
                25 / 81,
                40 / 81,
                25 / 81,
            ]
        else:
            raise ValueError("Not found npg = " + str(npg))
        return W


########################################################################
#             Definition of interpolation functions in 2D              #
########################################################################
#
#  We have two geometries:
#     Triangle:
#          . (0, 1)
#          |╲
#          | ╲
#          |  ╲
#          |   ╲
#          |____╲(1, 0)
#        (0,0)
#
#     Square:
#   (-1, 1)_________(1, 1)
#         |         |
#         |         |
#         |         |
#         |         |
#  (-1, -1)‾‾‾‾‾‾‾‾‾(1, -1)
#


class SFT3:
    """
    Shape functions for 3-node triangle
    Vertices:
        P1 = (0, 0)
        P2 = (1, 0)
        P3 = (0, 1)
    """

    @staticmethod
    def N(x):
        N1 = 1 - x[0] - x[1]
        N2 = x[0]
        N3 = x[1]
        return [N1, N2, N3]

    @staticmethod
    def dN(x):
        dN1_dx0 = -1
        dN1_dx1 = -1
        dN2_dx0 = 1
        dN2_dx1 = 0
        dN3_dx0 = 0
        dN3_dx1 = 1
        dN1 = [dN1_dx0, dN1_dx1]
        dN2 = [dN2_dx0, dN2_dx1]
        dN3 = [dN3_dx0, dN3_dx1]
        return [dN1, dN2, dN3]

    @staticmethod
    def P():
        P1 = (0, 0)
        P2 = (1, 0)
        P3 = (0, 1)
        return [P1, P2, P3]


class SFT6:
    """
    Shape functions for 3-node triangle
    Vertices:
        P1 = (0, 0)
        P2 = (0.5, 0)
        P3 = (1, 0)
        P4 = (0.5, 0.5)
        P5 = (0, 1)
        P6 = (0, 0.5)
    """

    @staticmethod
    def N(x):
        # Ni(Pj) = delta_ij
        N1 = (x[0] + x[1] - 1) * (2 * x[0] + 2 * x[1] - 1)
        N2 = 4 * x[0] * (1 - x[0] - x[1])
        N3 = x[0] * (2 * x[0] - 1)
        N4 = 4 * x[0] * x[1]
        N5 = x[1] * (2 * x[1] - 1)
        N6 = 4 * x[1] * (1 - x[0] - x[1])
        return [N1, N2, N3, N4, N5, N6]

    @staticmethod
    def dN(x):
        dN1 = [4 * (x[0] + x[1]) - 3, 4 * (x[0] + x[1]) - 3]
        dN2 = [4 * (1 - 2 * x[0] - x[1]), -4 * x[0]]
        dN3 = [4 * x[0] - 1, 0]
        dN4 = [4 * x[1], 4 * x[0]]
        dN5 = [0, 4 * x[1] - 1]
        dN6 = [-4 * x[1], 4 * (1 - x[0] - 2 * x[1])]
        return [dN1, dN2, dN3, dN4, dN5, dN6]

    @staticmethod
    def P():
        P1 = (0, 0)
        P2 = (0.5, 0)
        P3 = (1, 0)
        P4 = (0.5, 0.5)
        P5 = (0, 1)
        P6 = (0, 0.5)
        return [P1, P2, P3, P4, P5, P6]


class SFQ4:
    """
    Shape functions for 4-node quadrangle
    Vertices:
        P1 = (-1, -1)
        P2 = (1, -1)
        P3 = (1, 1)
        P4 = (-1, 1)
    """

    @staticmethod
    def N(x):
        N1 = (x[0] - 1) * (x[1] - 1) / 4
        N2 = (x[0] + 1) * (1 - x[1]) / 4
        N3 = (x[0] + 1) * (x[1] + 1) / 4
        N4 = (1 - x[0]) * (x[1] + 1) / 4
        return [N1, N2, N3, N4]

    @staticmethod
    def dN(x):
        dN1 = [(x[1] - 1) / 4, (x[0] - 1) / 4]
        dN2 = [(1 - x[1]) / 4, -(x[0] + 1) / 4]
        dN3 = [(x[1] + 1) / 4, (x[0] + 1) / 4]
        dN4 = [-(x[1] + 1) / 4, (1 - x[0]) / 4]
        return [dN1, dN2, dN3, dN4]

    @staticmethod
    def P():
        P1 = (-1, -1)
        P2 = (1, -1)
        P3 = (1, 1)
        P4 = (-1, 1)
        return [P1, P2, P3, P4]


class SFQ8:
    """
    Shape functions for 8-node quadrangle
    Vertices:
        P1 = (-1, -1)
        P2 = (0, -1)
        P3 = (1, -1)
        P4 = (1, 0)
        P5 = (1, 1)
        P6 = (0, 1)
        P7 = (-1, 1)
        P8 = (-1, 0)
    """

    @staticmethod
    def N(x):
        N1 = -(x[0] - 1) * (x[1] - 1) * (x[0] + x[1] + 1) / 4
        N2 = (x[0] - 1) * (x[0] + 1) * (x[1] - 1) / 2
        N3 = -(x[0] + 1) * (x[1] - 1) * (x[0] - x[1] - 1) / 4
        N4 = -(x[0] + 1) * (x[1] - 1) * (x[1] + 1) / 2
        N5 = (x[0] + 1) * (x[1] + 1) * (x[0] + x[1] - 1) / 4
        N6 = -(x[0] - 1) * (x[0] + 1) * (x[1] + 1) / 2
        N7 = (x[0] - 1) * (x[1] + 1) * (x[0] - x[1] + 1) / 4
        N8 = (x[0] - 1) * (x[1] - 1) * (x[1] + 1) / 2
        return [N1, N2, N3, N4, N5, N6, N7, N8]

    @staticmethod
    def dN(x):
        dN1 = [
            (2 * x[0] + x[1]) * (1 - x[1]) / 4,
            (1 - x[0]) * (x[0] + 2 * x[1]) / 4,
        ]
        dN2 = [x[0] * (x[1] - 1), (x[0] ** 2 - 1) / 2]
        dN3 = [
            (-2 * x[0] + x[1]) * (x[1] - 1) / 4,
            (-x[0] + 2 * x[1]) * (x[0] + 1) / 4,
        ]
        dN4 = [(1 - x[1] ** 2) / 2, -x[1] * (x[0] + 1)]
        dN5 = [
            (2 * x[0] + x[1]) * (x[1] + 1) / 4,
            (x[0] + 1) * (x[0] + 2 * x[1]) / 4,
        ]
        dN6 = [-x[0] * (x[1] + 1), (1 - x[0] ** 2) / 2]
        dN7 = [
            (2 * x[0] - x[1]) * (x[1] + 1) / 4,
            (x[0] - 1) * (x[0] - 2 * x[1]) / 4,
        ]
        dN8 = [(x[1] ** 2 - 1) / 2, x[1] * (x[0] - 1)]
        return [dN1, dN2, dN3, dN4, dN5, dN6, dN7, dN8]

    @staticmethod
    def P():
        P1 = (-1, -1)
        P2 = (0, -1)
        P3 = (1, -1)
        P4 = (1, 0)
        P5 = (1, 1)
        P6 = (0, 1)
        P7 = (-1, 1)
        P8 = (-1, 0)
        return [P1, P2, P3, P4, P5, P6, P7, P8]
