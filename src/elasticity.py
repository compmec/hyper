# -*- coding: utf-8 -*-
#
# Data structures for hyperelastic constitutive models
#
"""
#Potential, stress, stiffness for St Venant Kirchhoff / Neo-Hookean models.
WARNING:
in FiniteElement.update(), (S,M) = material.stress_stiffness(C)
--> stress_stiffness must be a function of right Cauchy Green tensor
===============================================================================
"""

import tensor
import numpy as np

class StVenantKirchhoffElasticity:
    """Data structure for (isotropic) St-Venant-Kirchhoff hyperelaticity models"""
    prop = dict()
    def __init__(self,E,nu):
        self.prop["YOUNG_MODULUS"] = E
        self.prop["POISSON_COEFFICIENT"] = nu
        self.prop["SHEAR_MODULUS"] = 0.5*E/(1.+nu)
        self.prop["BULK_MODULUS"] = E/3./(1.-2*nu)
        self.prop["1ST_LAME_CONSTANT"] = self.prop["BULK_MODULUS"]-2./3.*self.prop["SHEAR_MODULUS"]
        self.prop["2ND_LAME_CONSTANT"] = self.prop["SHEAR_MODULUS"]

    def potential(self,C):
        """Compute hyperelastic potential: phi = lambda/2 * tr(E)^2 - mu*(E:E)"""
        lam = self.prop["1ST_LAME_CONSTANT"]
        mu = self.prop["2ND_LAME_CONSTANT"]
        EL = 0.5*(C - tensor.I(len(C))) # Lagrangian strain E
        phi = lam/2.*(tensor.trace(EL))**2 + mu*np.tensordot(EL,EL,2)
        return phi

    def stress(self,C):
        """Compute 2nd Piola-Kirchhoff stress"""
        # =====================================================================
        # TODO: compute PK2, 2nd Piola-Kirchhoff stress tensor as a function of
        # C, right Cauchy Green tensor
        d = len(C)
        PK2 = tensor.tensor(d)
        lam = self.prop["1ST_LAME_CONSTANT"]
        mu = self.prop["2ND_LAME_CONSTANT"]
        E = tensor.tensor(d)
        E = 0.5*np.add(C,-tensor.I(d))
        PK2 = np.add(np.multiply(lam*tensor.trace(E),tensor.I(d)),np.multiply(2*mu,E))
        # =====================================================================
        return PK2

    def stiffness(self,C):
        """Compute material tangent M = 2*dS/dC """
        # =====================================================================
        # TODO: compute M, the material tangent as a function of
        # C, right Cauchy Green tensor
        d = len(C)
        lam = self.prop["1ST_LAME_CONSTANT"]
        mu = self.prop["2ND_LAME_CONSTANT"]
        M = tensor.tensor4(d)
        I = tensor.I(d)
        terme1 = tensor.outerProd4(I,I)
        terme1 = np.multiply(lam,terme1)
        terme2 = np.multiply(2*mu,tensor.IISym(d))
        M = np.add(terme1,terme2)
        # =====================================================================
        return M

    def stress_stiffness(self,C):
        """Compute 2nd Piola-Kirchhoff stress and material tangent at the same time"""
        # =====================================================================
        # TODO: compute PK2, 2nd Piola-Kirchhoff stress tensor as a function of
        # C, right Cauchy Green tensor
        PK2 = StVenantKirchhoffElasticity.stress(self,C)
        M = StVenantKirchhoffElasticity.stiffness(self,C)
        # =====================================================================
        return (PK2,M)

class NeoHookeanElasticity:
    """Data structure for (isotropic) NeoHookean hyperelaticity models"""
    prop = dict()
    def __init__(self,E,nu):
        self.prop["YOUNG_MODULUS"] = E
        self.prop["POISSON_COEFFICIENT"] = nu
        self.prop["SHEAR_MODULUS"] = 0.5*E/(1.+nu)
        self.prop["BULK_MODULUS"] = E/3./(1.-2*nu)
        self.prop["1ST_LAME_CONSTANT"] = self.prop["BULK_MODULUS"]-2./3.*self.prop["SHEAR_MODULUS"]
        self.prop["2ND_LAME_CONSTANT"] = self.prop["SHEAR_MODULUS"]

    def potential(self,C):
        """Compute hyperelastic potential: phi = mu/2 * (tr(C)-3) - mu*ln(J) + lam/2 *ln(J)^2"""
        lam = self.prop["1ST_LAME_CONSTANT"]
        mu = self.prop["2ND_LAME_CONSTANT"]
        J = np.sqrt(tensor.det(C)) # J = det(F) and det(C) = J^2
        phi = mu/2.*(tensor.trace(C)-3.) - mu*np.log(J) + lam/2.*(np.log(J))**2.
        return phi

    def stress(self,C):
        """Compute 2nd Piola-Kirchhoff stress"""
        # =====================================================================
        # TODO: compute PK2, 2nd Piola-Kirchhoff stress tensor as a function of
        # C, right Cauchy Green tensor
        d = len(C)
        PK2 = tensor.tensor(d)
        lnJ = np.log(np.sqrt(tensor.det(C)))
        lam = self.prop["1ST_LAME_CONSTANT"]
        mu = self.prop["2ND_LAME_CONSTANT"]
        invC = tensor.inv(C)
        terme1 = mu*np.add(tensor.I(d),-1.*invC)
        terme2 = lam * lnJ * invC
        PK2 = np.add(terme1,terme2)
        # =====================================================================
        return PK2

    def stiffness(self,C):
        """Compute material tangent M = 2*dS/dC """
        # =====================================================================
        # TODO: compute M, the material tangent as a function of
        # C, right Cauchy Green tensor
        d = len(C)
        M = tensor.tensor4(d)
        lnJ = np.log(np.sqrt(tensor.det(C)))
        lam = self.prop["1ST_LAME_CONSTANT"]
        mu = self.prop["2ND_LAME_CONSTANT"]
        invC = tensor.inv(C)
        terme1 = lam*tensor.outerProd4(invC,invC)
        dinvC = tensor.tensor4(d)
        for i in range(len(dinvC)):
            for j in range(len(dinvC)):
                for k in range(len(dinvC)):
                    for l in range(len(dinvC)):
                        dincC[i][j][k][l] = -0.5*(dinvC[i][k]*dinvC[j][l]+dinvC[i][l]*dinvC[j][k])


        terme2 = 2*(lam*lnJ-mu)*dinvC
        # =====================================================================
        return M

    def stress_stiffness(self,C):
        """Compute 2nd Piola-Kirchhoff stress and material tangent at the same time"""
        # =====================================================================
        # TODO: compute PK2, 2nd Piola-Kirchhoff stress tensor as a function of
        # C, right Cauchy Green tensor
        PK2 = NeoHookeanElasticity.stress(len(C))
        M = NeoHookeanElasticity.stiffness(len(C))
        # =====================================================================
        return (PK2,M)
