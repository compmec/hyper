# -*- coding: utf-8 -*-
#
# Data structures for finite element model
#


import numpy as np
import tensor
import geometry
from scipy.linalg import polar


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

        if t == 2:  # triangle T3 (1 int. pt.)
            npg = 3  # Number of integration points
            npN = 3  # Number of nodes in the element
            dim = 2
            IP = geometry.IPTri  # class to get Integration Info
            Shape = geometry.SFT3  # Function to get N(X) et dN(X)

        if t == 3:  # square Q4 (1 int. pt.)
            npg = 4  # Number of integration points
            npN = 4  # Number of nodes in the element
            dim = 2
            IP = geometry.IPQua  # class to get Integration Info
            Shape = geometry.SFQ4  # Function to get N(X) et dN(X)

        X = IP.getX(npg)
        W = IP.getW(npg)
        for i in range(npg):  # loop on integration points
            Xi = X[i]
            Wi = W[i]
            shapei = Shape.N(Xi)
            dshapei = Shape.dN(Xi)

            dshapei = np.array(dshapei)
            J = np.dot(np.transpose(xNod), dshapei)
            # for n in range(npN):  # loop on nodes
            #     xn = xNod[n][:dim]
            #     dshapein = dshapei[n]
            #     J += tensor.outerProd(xn, dshapein)
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
        d = self.getDim()  # space dimension
        npg = self.getNIntPts()
        for n in range(npg):
            self.F.append(tensor.tensor(d))
            self.hencky.append(tensor.tensor(d))
            self.E_GL.append(tensor.tensor(d))
            self.E_EA.append(tensor.tensor(d))
            self.PK1.append(tensor.tensor(d))
            self.sigma.append(tensor.tensor(d))
            self.K.append(tensor.tensor4(d))

    def getNpg(self):
        return np.shape(self.dshape)[0]

    def getNNodes(self):
        return np.shape(self.dshape)[1]

    def getDim(self):
        return np.shape(self.dshape)[2]

    def getNIntPts(self):
        return len(self.shape)

    def gradient(self, dShp, uNod):
        """
        Compute gradient of the displacement field
        """
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
        """
        update strain and stress tensors of the element from the displacement
        field 'uNod' with the material model 'mater'
        """
        # dshape = np.shape(self.dshape)
        # dshape = (npg, npN, dim)
        #   npg is the number of gauss points to integration
        #   npN is the number of nodes in each element
        #   dim is the dimention of the problem. If 2D -> dim = 2
        npg = self.getNpg()
        for i in range(npg):  # loop on integration points
            # compute CG stretch tensor and GL strain tensor
            # print("uNod = ")
            # print(uNod)
            # print("dshapei = ")
            # print(self.dshape[i])
            G = self.gradient(self.dshape[i], uNod)
            F = G + tensor.I(len(G))
            # print("G = ")
            # print(G)
            self.F[i] = F  # store deformation gradient at integration point i
            C = tensor.rightCauchyGreen(F)
            self.E_GL[i] = 0.5 * (C - tensor.I(len(C)))
            # compute spatial description
            # from scipy.linalg.polar() method with u=R and p=V,
            _, V = polar(F, side='left')
            # replace pure zeros by very low values to prevent "nan" in np.log(V)
            V[V < 1e-10] = 1e-15
            self.hencky[i] = np.log(V)  # ln(V) with F = V.R, "true" strain
            # print("F = ")
            # print(F)
            b = tensor.leftCauchyGreen(F)
            self.E_EA[i] = 0.5 * (tensor.I(len(b)) - tensor.inv(b))
            ###
            if (mater == 0):  # skip next lines: do not compute stress
                continue
            # compute PK2 stress tensor and material tangent operator M=2*dS/dC
            # print("type(mater) = " + str(type(mater)))
            # print("mater = " + str(mater))
            (PK2, M) = mater.stress_stiffness(C)
            # compute PK1 stress and lagrangian tangent operator K = dP/dF
            # print("PK2 = ")
            # print(PK2)
            # print("F = ")
            # print(F)
            self.PK1[i] = tensor.PK2toPK1(F, PK2)
            # print("PK1[i] = ")
            # print(self.PK1[i])
            self.K[i] = tensor.MaterialToLagrangian(F, PK2, M)
            # compute cauchy stress (spatial description)
            self.sigma[i] = tensor.PK1toCauchy(F, self.PK1[i])

    # def computeForces(self, fNod):
    def computeForces(self):
        """
        compute internal forces of the element
        """
        npg = self.getNpg()
        nNodes = self.getNNodes()
        dim = self.getDim()
        fNod = np.zeros((nNodes, dim))
        # dshape = (npg, npN, dim)
        #   npg is the number of gauss points to integration
        #   npN is the number of nodes in each element
        #   dim is the dimention of the problem. If 2D -> dim = 2
        for i in range(npg):
            PK1i = self.PK1[i]
            dshapei = self.dshape[i]
            weighti = self.weight[i]
            fNod += weighti * dshapei @ PK1i
        return fNod

    # def computeStiffness(self, KNod):
    def computeStiffness(self):
        """
        compute internal stiffness of the element
        """
        nNodes = self.getNNodes()
        dim = self.getDim()
        KNod = np.zeros((nNodes, dim, nNodes, dim))
        for i in range(npg):
            Ki = self.K[i]
            shapei = self.shape[i]
            KNod += weighti * Ki @ shapei
        return KNod


class FEModel:
    """
    Data structure for FE model
    """

    def __init__(self, theMesh, m=0):
        self.dim = theMesh.getDimension()  # mesh dimension
        self.elems = []  # list of FiniteElement instances
        self.connect = []  # table of connectivity (list of list)
        self.nnode = theMesh.getNNodes()  # Number of nodes in the mesh
        self.nelem = theMesh.getNElements()  # Number of elements in the mesh
        for n in range(self.nelem):
            elem = theMesh.getElement(n)
            if (elem.type < 1) or (elem.type > 3):
                continue
            self.connect.append(elem.nodes)
            xNod = []
            nNodes_elem = elem.getNNodes()
            for i in range(nNodes_elem):
                label_nodei = elem.nodes[i]
                nodei = theMesh.getNode(label_nodei)
                xnodei = nodei.getX()
                # Since we don't need z, we cut it out
                xnodei = xnodei[:self.dim]
                xNod.append(xnodei)
            xNod = np.array(xNod)
            new_FiniteElement = FiniteElement(elem.type, xNod)
            self.elems.append(new_FiniteElement)
        self.material = m  # constitutive model instance from elasticity module

    def getDim(self):
        """
        get the dimension of the mesh
        """
        return self.dim

    def getNNodes(self):
        """
        get number of nodes in the mesh
        """
        return self.nnode

    def getNElements(self):
        """
        get number of elements
        """
        return self.nelem

    def nElements(self):
        return self.getNElements()

    def assemble2(self, e, vNod):
        """
        Assemble nodal values for element e, 2-entry array
        """
        loc = self.connect[e]
        self.Tint[loc, :] += vNod[:, :]
        # raise Exception("assemble2: Not implemented")

    def assemble4(self, e, KNod):
        """
        Assemble nodal values for element e, 4-entry array
        """
        self.Kint[KNod, :, KNod, :] += K[:, :, :, :]
        # raise Exception("assemble4: Not implemented")

    def extract(self, U, e):
        """
        Extract nodal values for element e
        """
        if not isinstance(e, (int)):
            raise Exception(
                "The element 'e' in 'extract' function is not a integer")
        # print("index = " + str(index))
        connect_e = self.connect[e]
        # print("connect_e = " + str(connect_e))
        return [U[i] for i in connect_e]

    # def computeInternalForces(self, U, T):
    def computeInternalForces(self, U):
        """
        Compute generalized internal forces Tint
        """
        # T shape = (number of nodes in mesh, space dimension)
        # n, dim = T.shape
        nNodes = self.getNNodes()
        dim = self.getDim()
        nElements = self.getNElements()
        self.Tint = np.zeros((nNodes, dim))
        for e in range(nElements):
            uNod = self.extract(U, e)
            elem = self.elems[e]
            elem.update(uNod, self.material)
            vNod = elem.computeForces()
            self.assemble2(e, vNod)
        # print("Tint = ")
        # print(self.Tint)
        return self.Tint
        # raise Exception("computeInternalForces: Not implemented")

    # def computeResidual(self, U, R):
    def computeResidual(self, U):
        """
        Compute residual R = Tint - Text
        """
        nNodes = self.getNNodes()
        dim = self.getDim()
        R = np.zeros((nNodes, dim))
        R[:, :] = 0.0  # size = (number of nodes in mesh, space dimension)
        # R = self.computeInternalForces(U, R)
        R = self.computeInternalForces(U)  # R = Tint
        # R -= Text
        return R

    # def computeInternalStiffness(self, U, Kint):
    def computeInternalStiffness(self, U):
        """
        Compute generalized internal stiffness tangent Kint
        """
        # K shape = (number of nodes in mesh, space dimension,
        #            number of nodes in mesh, space dimension)
        nNodes = self.getNNodes()
        dim = self.getDim()
        self.Kint = np.zeros((nNodes, dim, nNodes, dim))
        for i, e in enumerate(self.elems):
            # uNod = self.extract(U, e)
            vNod = self.connect[i]
            K = e.computeStiffness()
            self.assemble4(e, vNod, K)
        print("Kint = ")
        print(self.Kint)
        return self.Kint
        #  size = (number of nodes in mesh, space dimension,number of nodes in mesh, space dimension)
        raise Exception("computeInternalStiffness: Not implemented")

    def computeTangent(self, U, K):
        """
        Compute tangent K = Kint - Kext = Kint: in the absence of displacement-
        depend external forces Kext = 0.
        """
        K[:, :, :, :] = 0.0
        #  size = (number of nodes in mesh, space dimension,number of nodes in mesh, space dimension)
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
        avg /= self.elems[n].getNIntPts()
        return avg

    def getStrainGreenLagrange(self, n):
        """
        Return average of GL strain at all integration points of element n
        """
        avg = tensor.tensor(self.dim)
        for tens in self.elems[n].E_GL:
            avg += tens
        avg /= self.elems[n].getNIntPts()
        return avg

    def getStressPK1(self, n):
        """
        Return average of PK1 stress at all integration points of element n
        """
        avg = tensor.tensor(self.dim)
        for tens in self.elems[n].PK1:
            avg += tens
        avg /= self.elems[n].getNIntPts()
        return avg

    def getStrainEulerAlmansi(self, n):
        """ Return average of EA strain at all integration points of element n"""
        avg = tensor.tensor(self.dim)
        for tens in self.elems[n].E_EA:
            avg += tens
        avg /= self.elems[n].getNIntPts()
        return avg

    def getStrainHencky(self, n):
        """
        Return average of hencky strain at all integration points of element n
        """
        avg = tensor.tensor(self.dim)
        for tens in self.elems[n].hencky:
            avg += tens
        avg /= self.elems[n].getNIntPts()
        return avg

    def getStressCauchy(self, n):
        """
        Return average of cauchy stress at all integration points of element n
        """
        avg = tensor.tensor(self.dim)
        for tens in self.elems[n].sigma:
            avg += tens
        avg /= self.elems[n].getNIntPts()
        return avg
