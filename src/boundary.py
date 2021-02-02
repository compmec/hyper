# -*- coding: utf-8 -*-

import numpy as np


class Boundary:
    """
    """

    def __init__(self, mesh):
        self.mesh = mesh
        dim = self.mesh.getDim()
        nNodes = self.mesh.getNNodes()
        nElem = self.mesh.getNElements()
        self.unitlength = 1e+9
        self.unitarea = 1e+9
        # we can say that unit is the minimum distance between two points
        #     unit = min dist(point_A - point_B)
        for i in range(nNodes):
            nodei = self.mesh.getNode(i)
            Xi = np.array(nodei.getX())
            for j in range(i + 1, nNodes):
                nodej = self.mesh.getNode(j)
                Xj = np.array(nodej.getX())
                dX = Xi - Xj
                if np.linalg.norm(dX) < self.unitlength:
                    self.unitlength = np.linalg.norm(dX)

        for i in range(nElem):
            elem = self.mesh.getElement(i)
            if elem.getType() == 2:  # Triangle
                X0 = np.array(self.mesh.getNode(elem.getNode(0)).getX())
                X1 = np.array(self.mesh.getNode(elem.getNode(1)).getX())
                X2 = np.array(self.mesh.getNode(elem.getNode(2)).getX())
                dX1 = X1 - X0
                dX2 = X2 - X0
                area = np.linalg.norm(np.cross(dX1, dX2))
                if area < self.unitarea:
                    self.unitarea = area
            elif elem.getType() == 3:  # Quadrangle
                pass

        self.PRECISION = self.unitlength / 100
        self.bcDofs = []
        self.freeDofs = []
        for i in range(nNodes):
            for j in range(dim):
                self.freeDofs.append(i * dim + j)

    def get_points(self, fs=[], gs=[]):
        """
        This function gets all the points that satisfy the vectorial function.
        For exemple, if the function is f(x, y) = x^2 + y^2 - R^2
        so, this function will return all the points such that
            f(x, y) == 0
        In this case, all the points near to a circle of radius R
        As we don't have always exact values, but approximate values, to implement the function we say that:
            abs( f(x, y) ) < PRECISION
        We can say that there is a group of functions inside gs.
        For exemple, if we want all the points in the frist quadrant, we have to say that
            g1(x, y) = x
            g2(x, y) = y
        And the restriction is that g >= 0
        To implement, the points will satisfy g > -PRECISION
        """
        nNodes = self.mesh.getNNodes()
        points = []
        for i in range(nNodes):
            nodei = self.mesh.getNode(i)
            Xi = nodei.getX()
            possible = True
            for f in fs:
                if np.abs(f(Xi)) > self.PRECISION:
                    possible = False
                    break
            if not possible:
                continue
            for g in gs:
                if g(Xi) < -self.PRECISION:
                    possible = False
                    break
            if not possible:
                continue
            points.append(i)
        return points

    def add_restriction(self, points, directions):
        dim = self.mesh.getDim()
        dofs = []
        if isinstance(directions, str):
            directions = [directions]
        for i, p in enumerate(points):
            if "x" in directions:
                dofs.append(p * dim)
            if "y" in directions:
                dofs.append(p * dim + 1)
            if "z" in directions:
                dofs.append(p * dim + 2)

        for d in dofs:
            if d not in self.bcDofs:
                self.bcDofs.append(d)
                self.freeDofs.remove(d)

        return dofs

    def cut(self, element, axis=0):
        if not isinstance(element, np.ndarray):
            raise Exception("Expected a numpy array to cut")
        if len(self.freeDofs) + len(self.bcDofs) != element.shape[axis]:
            msg = "Size of numpy array is not compatible.\n"
            msg += "bcDofs + freeDofs != shape --- %d + %d = %d != %d" % (len(self.bcDofs), len(
                self.freeDofs), len(self.freeDofs) + len(self.bcDofs), element.shape[0])
            raise Exception(msg)
        K1 = np.delete(element, self.freeDofs, axis=axis)
        K2 = np.delete(element, self.bcDofs, axis=axis)
        if len(element.shape) == axis + 1:
            return [K1, K2]
        else:
            return [self.cut(K1, axis + 1), self.cut(K2, axis + 1)]

    def getBcDofs(self):
        return self.bcDofs

    def getFreeDofs(self):
        return self.freeDofs

    def getUnitLength(self):
        return self.unitlength

    def getUnitArea(self):
        return self.unitarea
