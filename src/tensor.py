# -*- coding: utf-8 -*-
#
#        @author: Clément Fuchs
#                 Carlos Adir Ely Murussi Leite
#
#
# Set of methods for manipulating 2nd-order and 4th-order tensors
#

import numpy as np
import numpy.linalg as la


###############################################################################
#                                   Vectors                                   #
###############################################################################

def vector(d=2):
    """
    Constructor of a vector object (dimension d)
    """
    return np.zeros(d)


###############################################################################
#                             2nd Order Tensors                               #
###############################################################################

def tensor(d=2):
    """
    Constructor of 2nd-order tensor (dimension d)
    """
    return np.zeros((d, d))


def I(d=2):
    """
    Identity second-order tensor
    """
    return np.eye(d)


def det(A):
    """
    Compute determinant of a matrix/tensor
    """
    if 1e-13 < np.abs(la.det(A)) < 1e-10:
        print("A = ")
        print(A)
        print("detA = ")
        print(la.det(A))
        raise Exception("Not good")
    return la.det(A)


def inv(A):
    """
    Compute inverse of a matrix/tensor
    """
    return la.inv(A)


def trace(A):
    """
    Compute trace of a matrix/tensor
    """
    return np.trace(A)


def outerProd(a, b):
    """
    Compute outer product of two vectors
        Mij = ai * bj
    """
    return np.tensordot(a, b, axes=0)


def rightCauchyGreen(F):
    """
    Compute right Cauchy-Green tensor from deformation gradient
    """
    return np.dot(np.transpose(F), F)


def leftCauchyGreen(F):
    """
    Compute left Cauchy-Green tensor from deformation gradient
    """
    return np.dot(F, np.transpose(F))


def PK2toPK1(F, S):
    """
    Compute Piola stress tensor from second Piola-Kirchhoff stress
    """
    return np.dot(F, S)


def PK1toCauchy(F, P):
    """
    Compute Cauchy stress tensor from first Piola-Kirchhoff stress
    """
    return np.dot(P, np.transpose(F))

###############################################################################
#                             4th Order Tensors                               #
###############################################################################


def tensor4(d=2):
    """
    Constructor of 4th-order tensor (dimension d)
    """
    return np.zeros((d, d, d, d))


def II(d=2):
    """
    Identity fourth-order tensor
        II1_ijkl = delta_ik * delta_jl
    """
    II1 = tensor4(d)
    for i in range(d):
        for j in range(d):
            II1[i, j, i, j] = 1
    return II1


def JJ(d=2):
    """
    Transpose fourth-order tensor
        JJ1_ijkl = delta_il * delta_jk
    If we have the tensor A, so
        A^T = JJ : A
    """
    JJ1 = tensor4(d)
    for i in range(d):
        for j in range(d):
            JJ1[i, j, j, i] = 1
    return JJ1


def IISym(d=2):
    """
    Symmetrical identity fourth - order tensor:
        IISym_ijkl = 1 / 2 * (delta_ik delta_jl + delta_il delta_jk)
    If we have the 2nd order tensor A:
        IISym : A = (A + A^T)/2
    """
    return (II(d) + JJ(d)) / 2


def KK(d=2):
    """
    Spherical operator
    """
    S = tensor4(d)
    for i in range(d):
        for k in range(d):
            S[i, i, k, k] = 1
    return S


def outerProd4(a, b):
    """
    Compute outer product of two tensors
        M_ijkl = a_ij * b_kl
    """
    return np.tensordot(a, b, axes=0)


def MaterialToLagrangian(F, S, M):
    """
    Compute Lagrangian tangent operator from material tensor and stress
    """
    # raise Exception("MaterialToLagrangian: Not implemented")
    return "not implemented"
