# -*- coding: utf-8 -*-

import tensor
import numpy as np


class Elasticity:
    """
    Data structure for (isotropic) hyperelaticity models
    """
    prop = dict()

    def __init__(self, E, nu):
        self.prop["YOUNG_MODULUS"] = E
        self.prop["POISSON_COEFFICIENT"] = nu
        self.prop["SHEAR_MODULUS"] = 0.5 * E / (1. + nu)
        self.prop["BULK_MODULUS"] = E / 3. / (1. - 2 * nu)
        self.prop["1ST_lamdaE_CONSTANT"] = self.prop["BULK_MODULUS"] - \
            2. / 3. * self.prop["SHEAR_MODULUS"]
        self.prop["2ND_lamdaE_CONSTANT"] = self.prop["SHEAR_MODULUS"]

    def stress_stiffness(self, C):
        """
        Compute 2nd Piola-Kirchhoff stress and material tangent at the same time
        """
        PK2 = self.stress(C)
        M = self.stiffness(C)
        return (PK2, M)


class StVenantKirchhoffElasticity(Elasticity):

    def __init__(self, E, nu):
        super(StVenantKirchhoffElasticity, self).__init__(E, nu)

    def potential(self, C):
        """
        Compute hyperelastic potential: phi = lamdabda/2 * tr(E)^2 - mu*(E:E)
        """
        lamda = self.prop["1ST_lamdaE_CONSTANT"]
        mu = self.prop["2ND_lamdaE_CONSTANT"]
        EL = 0.5 * (C - tensor.I(len(C)))  # Lagrangian strain E
        phi = lamda / 2. * (tensor.trace(EL))**2 + mu * np.tensordot(EL, EL, 2)
        return phi

    def stress(self, C):
        """
        Compute 2nd Piola-Kirchhoff stress
        """
        d = len(C)
        PK2 = tensor.tensor(d)
        lamda = self.prop["1ST_lamdaE_CONSTANT"]
        mu = self.prop["2ND_lamdaE_CONSTANT"]
        E = tensor.tensor(d)
        E = 0.5 * (C - tensor.I(d))
        II = tensor.I(d)
        PK2 = lamda * tensor.trace(E) * II + 2 * mu * E
        return PK2

    def stiffness(self, C):
        """
        Compute material tangent M = 2*dS/dC
        """
        d = len(C)
        lamda = self.prop["1ST_lamdaE_CONSTANT"]
        mu = self.prop["2ND_lamdaE_CONSTANT"]
        I = tensor.I(d)
        IxI = tensor.outerProd4(I, I)
        IIsym = tensor.IISym(d)
        M = lamda * IxI + 2 * mu * IIsym
        return M


class NeoHookeanElasticity(Elasticity):

    def __init__(self, E, nu):
        super(NeoHookeanElasticity, self).__init__(E, nu)

    def potential(self, C):
        """
        Compute hyperelastic potential: phi = mu/2 * (tr(C)-3) - mu*ln(J) + lamda/2 *ln(J)^2
        """
        lamda = self.prop["1ST_lamdaE_CONSTANT"]
        mu = self.prop["2ND_lamdaE_CONSTANT"]
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
        d = len(C)
        PK2 = tensor.tensor(d)
        detC = tensor.det(C)
        detF = np.sqrt(detC)
        lnJ = np.log(detF)
        lamda = self.prop["1ST_lamdaE_CONSTANT"]
        mu = self.prop["2ND_lamdaE_CONSTANT"]
        invC = tensor.inv(C)
        I = tensor.I(d)
        PK2 = mu * (I - invC) + lamda * lnJ * invC
        return PK2

    def stiffness(self, C):
        """
        Compute material tangent M = 2*dS/dC
        """
        d = len(C)
        lnJ = np.log(np.sqrt(tensor.det(C)))
        lamda = self.prop["1ST_lamdaE_CONSTANT"]
        mu = self.prop["2ND_lamdaE_CONSTANT"]
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
