# -*- coding: utf-8 -*-
#
# Data structures for unstructured mesh
#
"""
MNMNL Homework
Module 2: nodal field on a mesh
complete code below where indicated by

----------------------------

Here we have all the definitions of the base elements:
* Nodes
* Elements
* Mesh
And we are going to use it futher

"""


class Node:
    """Basic data structure for nodes"""

    def __init__(self, i, x, y, z):
        """Create node of label i and coordinates (x,y,z)"""
        self.id = i
        self.X = (x, y, z)

    def getID(self):
        return self.id

    def getX(self, i=None):
        if i is None:
            return self.X
        else:
            return self.X[i]


class Element:
    """Basic data structure for elements"""

    def __init__(self, i, t, n):
        """Create element of label i, type t, with list of nodes n"""
        self.id = i
        self.type = t    # 2:Triangle, 3:Quadrangle
        self.nodes = []  # list of nodes labels (int)
        self.nodes.extend(n)

    def getID(self):
        return self.id

    def getType(self):
        return self.type

    def nNodes(self):
        return len(self.nodes)

    def getNode(self, i):
        return self.nodes[i]


class MeshData:
    """Class containing basic data structure for unstructured mesh"""

    def __init__(self, d=2):
        """Create mesh of dimension d"""
        self.dim = d  # dimension of the mesh as an integer, 2 by default
        self.nodes = []  # list of all the nodes in the mesh as Node instances
        self.elems = []  # list of all the elements in the mesh as Element instances

    def getDimension(self):
        # =====================================================================
        # Returns the dimention of the Mesh
        return self.dim

    def addNode(self, i, x, y=0., z=0.):
        self.nodes.append(Node(i, x, y, z))

    def getNode(self, i):
        # =====================================================================
        return self.nodes[i]
        # =====================================================================

    def nNodes(self):
        # =====================================================================
        # TODO: return number of nodes in the mesh
        return len(self.nodes)
        # =====================================================================

    def addElement(self, i, t, n):
        # =====================================================================
        # TODO: add element with label i, of type t, with list of nodes labels n
        new_element = Element(i, t, n)
        return self.elems.append(new_element)
        # =====================================================================

    def getElement(self, i):
        # =====================================================================
        # TODO: return element with label i
        return self.elems[i]
        # =====================================================================

    def nElements(self):
        # =====================================================================
        # TODO: return number of elements in the mesh
        return len(self.elems)
        # =====================================================================


if __name__ == "__main__":
    print("You should use it as a libray, not main file")
