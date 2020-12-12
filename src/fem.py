# -*- coding: utf-8 -*-
#
# Data structures for finite element model
#


import numpy as np
import geometrie
from scipy.linalg import polar

import tensor


########################################################################
#          Definition of elements for Finite Elements Method           #
########################################################################

class FiniteElement:
    """
    Data structure for finite element
    """

    def __init__(self, t, xNod):
        self.type = t
        self.shape = []
        self.dshape = []
        self.weight = []
        #
        # --- select element type and relevant shape functions
        #
        if t == 2:  # triangle T3 (1 int. pt.)
            npg = 3  # Number of integration points
            npN = 3  # Number of nodes in the element
            IP = geometrie.IPTri  # class to get Integration Info
            Shape = geometrie.SFT3  # Function to get N(X) et dN(X)

        if t == 3:  # square Q4 (1 int. pt.)
            npg = 4  # Number of integration points
            npN = 4  # Number of nodes in the element
            IP = geometrie.IPQua  # class to get Integration Info
            Shape = geometrie.SFQ4  # Function to get N(X) et dN(X)

        for i in range(npg):  # loop on integration points
            Xi = IP.X[i]
            Wi = IP.W[i]
            shapei = Shape.N(Xi)
            dshapei = Shape.dN(Xi)

            dshapei = np.array(dshapei)
            J = tensor.tensor(2)  # dX/dKsi
            for n in range(npN):  # loop on nodes
                J += tensor.outerProd(xNod[n][0:2], dshapei[n])
            detJ = tensor.det(J)
            Jinv = tensor.inv(J)
            dshapei = np.dot(dshapei, Jinv)
            weighti = Wi * detJ
            self.shape.append(shapei)
            self.dshape.append(dshapei)
            self.weight.append(weighti)

        #
        # --- initialise mechanical tensors at each integration point
        #
        self.F = []        # deformation gradient F = Grad(u) + I
        self.hencky = []   # hencky strain ln(V) with F = V.R
        self.E_GL = []     # lagrangian strain E = 1/2 (C - I)
        self.E_EA = []     # E_EA strain e = 1/2 (I - b^-1)
        self.PK1 = []      # piola kirchoff I : P = F*S
        self.sigma = []    # cauchy stress : \sigma
        self.K = []        # lagrangian tangent operator dP/dF
        d = self.getdim()  # space dimension
        for n in range(self.getnIntPts()):
            self.F.append(tensor.tensor(d))
            self.hencky.append(tensor.tensor(d))
            self.E_GL.append(tensor.tensor(d))
            self.E_EA.append(tensor.tensor(d))
            self.PK1.append(tensor.tensor(d))
            self.sigma.append(tensor.tensor(d))
            self.K.append(tensor.tensor4(d))

    def getdim(self,):
        return np.shape(self.dshape)[2]

    def getnIntPts(self,):
        return len(self.shape)

    def getnNodes(self,):
        return np.shape(self.dshape)[1]

    def gradient(self, dShp, uNod):
        """Compute gradient of the displacement field"""
        # =====================================================================
        # TODO: compute the gradient of the displacement field from the local
        # nodal field 'uNod' and the shape functions 'dShp' (20pts)
        # ...
        uNod = np.array(uNod)
        dShp = np.transpose(dShp)
        grad = np.dot(dShp, uNod)
        return grad
        # =====================================================================

    def update(self, uNod, mater):
        """ update strain and stress tensors of the element from the displacement
        field 'uNod' with the material model 'mater'"""
        for i in range(np.shape(self.dshape)[0]):  # loop on integration points
            # compute CG stretch tensor and GL strain tensor
            G = self.gradient(self.dshape[i], uNod)
            F = G + tensor.I(len(G))
            self.F[i] = F  # store deformation gradient at integration point i
            C = tensor.rightCauchyGreen(F)
            self.E_GL[i] = 0.5 * (C - tensor.I(len(C)))
            # compute spatial description
            # from scipy.linalg.polar() method with u=R and p=V,
            _, V = polar(F, side='left')
            # replace pure zeros by very low values to prevent "nan" in np.log(V)
            V[V < 1e-10] = 1e-15
            self.hencky[i] = np.log(V)  # ln(V) with F = V.R, "true" strain
            b = tensor.leftCauchyGreen(F)
            self.E_EA[i] = 0.5 * (tensor.I(len(b)) - tensor.inv(b))
            ###
            if (mater == 0):  # skip next lines: do not compute stress
                continue
            # compute PK2 stress tensor and material tangent operator M=2*dS/dC
            (PK2, M) = mater.stress_stiffness(C)
            # compute PK1 stress and lagrangian tangent operator K = dP/dF
            self.PK1[i] = tensor.PK2toPK1(F, PK2)
            self.K[i] = tensor.MaterialToLagrangian(F, PK2, M)
            # compute cauchy stress (spatial description)
            self.sigma[i] = tensor.PK1toCauchy(F, self.PK1[i])

    def computeForces(self, fNod):
        """ compute internal forces of the element"""
        fNod[:,
             :] = 0.0  # size = (number of nodes per element, space dimension)
        raise Exception("computeForces: Not implemented")

    def computeStiffness(self, KNod):
        """ compute internal stiffness of the element"""
        KNod[:, :, :, :] = 0.0  # shape : (nNodes, dim, nNodes, dim)
        raise Exception("computeStiffness: Not implemented")


class FEModel:
    """
    Data structure for FE model
    """

    def __init__(self, theMesh, m=0):
        self.dim = theMesh.getDimension()  # mesh dimension
        self.elems = []  # list of FiniteElement instances
        self.connect = []  # table of connectivity (list of list)
        for n in range(theMesh.nElements()):
            if (theMesh.elems[n].type < 1) or (theMesh.elems[n].type > 3):
                continue
            self.connect.append(theMesh.elems[n].nodes)
            xNod = []
            for i in range(theMesh.elems[n].nNodes()):
                xNod.append(theMesh.nodes[theMesh.elems[n].nodes[i]].X)
            new_FiniteElement = FiniteElement(theMesh.elems[n].type, xNod)
            self.elems.append(new_FiniteElement)
        self.material = m  # constitutive model instance from elasticity module

    def nElements(self):
        """
        Get number of elements
        """
        return len(self.elems)

    def assemble2(self, e, vNod, V):
        """
        Assemble nodal values for element e, 2-entry array
        """
        raise Exception("assemble2: Not implemented")

    def assemble4(self, e, KNod, K):
        """
        Assemble nodal values for element e, 4-entry array
        """
        raise Exception("assemble4: Not implemented")

    def extract(self, U, e):
        """
        Extract nodal values for element e
        """
        connect_e = self.connect[e]
        return [U[i] for i in connect_e]

    def computeInternalForces(self, U, T):
        """
        Compute generalized internal forces Tint
        """
        T[:, :] = 0.0  # size = (number of nodes in mesh, space dimension)
        raise Exception("computeInternalForces: Not implemented")

    def computeResidual(self, U, R):
        """
        Compute residual R = Tint - Text
        """
        R[:, :] = 0.0  # size = (number of nodes in mesh, space dimension)
        self.computeInternalForces(U, R)  # R = Tint
        return R

    def computeInternalStiffness(self, U, Kint):
        """
        Compute generalized internal stiffness tangent Kint
        """
        Kint[:, :, :,
             :] = 0.0  # size = (number of nodes in mesh, space dimension,number of nodes in mesh, space dimension)
        raise Exception("computeInternalStiffness: Not implemented")

    def computeTangent(self, U, K):
        """
        Compute tangent K = Kint - Kext = Kint: in the absence of displacement-
        depend external forces Kext = 0.
        """
        K[:, :, :,
            :] = 0.0  # size = (number of nodes in mesh, space dimension,number of nodes in mesh, space dimension)
        self.computeInternalStiffness(U, K)  # K = Kint
        return K

    def computeStrain(self, U):
        """
        Compute strains on all elements
        """
        # loop on elements
        for n in range(len(self.elems)):
            # get nodal displacements for element
            uNod = self.extract(U, n)
            # compute gradients (and strains)
            # print("uNod = " + str(uNod))
            self.elems[n].update(uNod, mater=0)

    def computeStress(self, U):
        """
        Compute stresses on all elements
        """
        # loop on elements
        for n in range(len(self.elems)):
            # get nodal displacements for element
            uNod = self.extract(U, n)
            # compute stresses
            self.elems[n].update(uNod, self.material)

    def getDeformationGradient(self, n):
        """
        Return average of deformation gradient at all integration points of element n
        """
        avg = tensor.tensor(self.dim)
        for tens in self.elems[n].F:
            avg += tens
        avg /= self.elems[n].getnIntPts()
        return avg

    def getStrainGreenLagrange(self, n):
        """
        Return average of GL strain at all integration points of element n
        """
        avg = tensor.tensor(self.dim)
        for tens in self.elems[n].E_GL:
            avg += tens
        avg /= self.elems[n].getnIntPts()
        return avg

    def getStressPK1(self, n):
        """
        Return average of PK1 stress at all integration points of element n
        """
        avg = tensor.tensor(self.dim)
        for tens in self.elems[n].PK1:
            avg += tens
        avg /= self.elems[n].getnIntPts()
        return avg

    def getStrainEulerAlmansi(self, n):
        """ Return average of EA strain at all integration points of element n"""
        avg = tensor.tensor(self.dim)
        for tens in self.elems[n].E_EA:
            avg += tens
        avg /= self.elems[n].getnIntPts()
        return avg

    def getStrainHencky(self, n):
        """
        Return average of hencky strain at all integration points of element n
        """
        avg = tensor.tensor(self.dim)
        for tens in self.elems[n].hencky:
            avg += tens
        avg /= self.elems[n].getnIntPts()
        return avg

    def getStressCauchy(self, n):
        """
        Return average of cauchy stress at all integration points of element n
        """
        avg = tensor.tensor(self.dim)
        for tens in self.elems[n].sigma:
            avg += tens
        avg /= self.elems[n].getnIntPts()
        return avg
