# -*- coding: utf-8 -*-

from abc import ABC, abstractmethod  # For abstract classes
from . import tensor
import numpy as np


class Elasticity(ABC):
    """
    Data structure for (isotropic) hyperelaticity models
    1ST_lamdaE_CONSTANT is the frist lamé parameter lambda
    2ST_lamdaE_CONSTANT is the second lamé parameter mu
    """
    prop = dict()

    def __init__(self, E, nu):
        G = E / (2 * (1. + nu))
        K = E / (3 * (1 - 2 * nu))
        lamda = K - (2 / 3) * G
        mu = G
        self.prop["YOUNG_MODULUS"] = E
        self.prop["POISSON_COEFFICIENT"] = nu
        self.prop["SHEAR_MODULUS"] = G
        self.prop["BULK_MODULUS"] = K
        self.prop["1ST_LAME_CONSTANT"] = lamda
        self.prop["2ND_LAME_CONSTANT"] = mu

    def stress_stiffness(self, C):
        """
        Compute 2nd Piola-Kirchhoff stress and material tangent at the same time
        """
        PK2 = self.stress(C)
        M = self.stiffness(C)
        return (PK2, M)

    @abstractmethod
    def potential(self, C):
        raise Exception("Elasticity - Potential Function not implemented!")

    @abstractmethod
    def stress(self, C):
        raise Exception("Elasticity - Stress Function not implemented!")

    @abstractmethod
    def stiffness(self, C):
        raise Exception("Elasticity - Stiffness Function not implemented!")

    @property
    def G(self):
        return self.prop["SHEAR_MODULUS"]

    @property
    def E(self):
        return self.prop["YOUNG_MODULUS"]

    @property
    def K(self):
        return self.prop["BULK_MODULUS"]

    @property
    def nu(self):
        return self.prop["POISSON_COEFFICIENT"]

    @property
    def LAME1(self):
        return self.prop["1ST_LAME_CONSTANT"]

    @property
    def LAME2(self):
        return self.prop["2ND_LAME_CONSTANT"]

    def get1LAME(self):
        return self.LAME1

    def get2LAME(self):
        return self.LAME2


class StVenantKirchhoffElasticity(Elasticity):

    def __init__(self, E, nu):
        super(StVenantKirchhoffElasticity, self).__init__(E, nu)

    def potential(self, C):
        """
        Compute hyperelastic potential: phi = lamdabda/2 * tr(E)^2 - mu*(E:E)
        """
        lamda = self.get1LAME()
        mu = self.get2LAME()
        EL = 0.5 * (C - tensor.I(len(C)))  # Lagrangian strain E
        phi = lamda / 2. * (tensor.trace(EL))**2 + mu * np.tensordot(EL, EL, 2)
        return phi

    def stress(self, C):
        """
        Compute 2nd Piola-Kirchhoff stress
        PK2 = 2*dphi/dC
            = lambda * tr(E) * I + 2 * mu * E
        """
        dim = len(C)
        PK2 = tensor.tensor(dim)
        lamda = self.get1LAME()
        mu = self.get2LAME()
        E = tensor.tensor(dim)
        E = 0.5 * (C - tensor.I(dim))
        II = tensor.I(dim)
        PK2 = lamda * tensor.trace(E) * II + 2 * mu * E
        return PK2

    def stiffness(self, C):
        """
        Compute material tangent M = 2*dS/dC
        """
        dim = len(C)
        lamda = self.get1LAME()
        mu = self.get2LAME()
        I = tensor.I(dim)
        IxI = tensor.outerProd4(I, I)
        IIsym = tensor.IISym(dim)
        M = lamda * IxI + 2 * mu * IIsym
        return M


class NeoHookeanElasticity(Elasticity):

    def __init__(self, E, nu):
        # self.VERBOSE = True
        super(NeoHookeanElasticity, self).__init__(E, nu)

    def potential(self, C):
        """
        Compute hyperelastic potential: phi = mu/2 * (tr(C)-3) - mu*ln(J) + lamda/2 *ln(J)^2
        """
        lamda = self.get1LAME()
        mu = self.get2LAME()
        dim = len(C)
        if dim == 2:
            K = np.copy(C)
            C = np.zeros((3, 3))
            C[:2, :2] = K[:, :]
        J = np.sqrt(tensor.det(C))  # J = det(F) and det(C) = J^2
        part1 = (mu / 2) * (tensor.trace(C) - 3.)
        part2 = mu * np.log(J)
        part3 = (lamda / 2) * (np.log(J))**2
        phi = part1 - part2 + part3
        return phi

    def stress(self, C):
        """
        Compute 2nd Piola-Kirchhoff stress
        """

        dim = len(C)
        if dim == 2:
            K = np.copy(C)
            C = np.eye(3)
            C[:2, :2] = K[:, :]
        PK2 = tensor.tensor(dim)
        detC = tensor.det(C)
        detF = np.sqrt(detC)
        lnJ = np.log(detF)
        lamda = self.get1LAME()
        mu = self.get2LAME()
        invC = tensor.inv(C)
        I = tensor.I(3)
        PK2 = mu * (I - invC) + lamda * lnJ * invC
        if dim == 2:
            return PK2[:2, :2]
        return PK2

    def stiffness(self, C):
        """
        Compute material tangent M = 2*dS/dC
        """
        d = len(C)
        lnJ = np.log(np.sqrt(tensor.det(C)))
        lamda = self.get1LAME()
        mu = self.get2LAME()
        invC = tensor.inv(C)
        invCC = tensor.outerProd4(invC, invC)
        terme1 = lamda * invCC
        dinvC = tensor.tensor4(d)
        for i in range(d):
            for j in range(d):
                for k in range(d):
                    for l in range(d):
                        part1 = invC[i, k] * invC[j, l]
                        part2 = invC[i, l] * invC[j, k]
                        dinvC[i, j, k, l] = -(part1 + part2) / 2

        M = lamda * invCC + 2 * (lamda * lnJ - mu) * dinvC

        return M
