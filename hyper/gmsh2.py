# -*- coding: utf-8 -*-
#
# Functions to write Gmsh output files
#
import numpy as np
import hyper.mesh as msh

#
# Read mesh data from a .msh file in Gmsh V.2 or V.4 format
# INPUT: - inp : input stream
#        - d : spatial dimension required for mesh data object (default = 2)
#              NOTE: all Gmsh mesh files are 3d
# OUTPUT: a MeshData object (as defined in mesh.py)
#


def gmshInput_mesh(inp, d=2):
    """Read mesh in Gmsh format"""
    lines = inp.readlines()
    mesh = msh.MeshData()
    iTab = dict()  # table matching node labels to internal indices
    #
    # --- determine format (v.2 or v.4 more recent)
    #
    words = lines[1].split()  # $MeshFormat information
    fmt = int(float(words[0]))  # version format 2 or 4
    fbin = int(words[1])  # whether binary file or not
    llen = int(words[2])  # data-size, unused
    # loop on lines
    l = 0
    #
    # --- version 2
    #
    if fmt == 2:
        while (l < len(lines)):
            if lines[l].startswith("$Nodes"):
                l += 1
                nNodes = int(lines[l].strip())
                i0 = len(iTab)
                for n in range(nNodes):
                    l += 1
                    words = lines[l].split()
                    lbl = int(words[0])
                    iTab[lbl] = i0 + n
                    if d == 1:
                        mesh.addNode(lbl, float(words[1]))
                    elif d == 2:
                        mesh.addNode(lbl, float(words[1]), float(words[2]))
                    elif d == 3:
                        mesh.addNode(lbl, float(words[1]), float(
                            words[2]), float(words[3]))
            elif lines[l].startswith("$Elements"):
                l += 1
                nElems = int(lines[l].strip())
                for n in range(nElems):
                    l += 1
                    words = lines[l].split()
                    lbl = int(words[0])
                    typ = int(words[1])
                    ntg = int(words[2])
                    iNodes = []
                    if (d == 1 and typ == 1):  # 2-node edge
                        iNodes.append(iTab[int(words[ntg + 3])])
                        iNodes.append(iTab[int(words[ntg + 4])])
                        mesh.addElement(lbl, typ, iNodes)
                    elif (d == 2 and typ == 2):  # 3-node triangle
                        iNodes.append(iTab[int(words[ntg + 3])])
                        iNodes.append(iTab[int(words[ntg + 4])])
                        iNodes.append(iTab[int(words[ntg + 5])])
                        mesh.addElement(lbl, typ, iNodes)
                    elif (d == 2 and typ == 3):  # 4-node quadrangle
                        iNodes.append(iTab[int(words[ntg + 3])])
                        iNodes.append(iTab[int(words[ntg + 4])])
                        iNodes.append(iTab[int(words[ntg + 5])])
                        iNodes.append(iTab[int(words[ntg + 6])])
                        mesh.addElement(lbl, typ, iNodes)
            else:
                l += 1  # ignore line
    #
    # --- version 4
    #
    elif fmt == 4:
        while (l < len(lines)):
            if lines[l].startswith("$Nodes"):
                l += 1
                words = lines[l].split()
                nEntityBlocks = int(words[0])
                nNodes = int(words[1])
                i0 = len(iTab)
                numNode = -1
                for n in range(nEntityBlocks):
                    l += 1
                    words = lines[l].split()
                    tagE = int(words[0])
                    dimE = int(words[1])
                    param = int(words[2])  # if parametric, 1, else 0
                    if param != 0:
                        # != instead of <> in Python 3
                        raise ValueError(
                            'parametric nodes are not handled for now')
                    numNodes = int(words[3])
                    for n in range(numNodes):
                        l += 1
                        words = lines[l].split()
                        lbl = int(words[0])
                        numNode += 1
                        iTab[lbl] = i0 + numNode
                        if d == 1:
                            mesh.addNode(lbl, float(words[1]))
                        elif d == 2:
                            mesh.addNode(lbl, float(words[1]), float(words[2]))
                        elif d == 3:
                            mesh.addNode(lbl, float(words[1]), float(
                                words[2]), float(words[3]))
            elif lines[l].startswith("$Elements"):
                l += 1
                words = lines[l].split()
                nEntityBlocks = int(words[0])
                nElems = int(words[1])
                for n in range(nEntityBlocks):
                    l += 1
                    words = lines[l].split()
                    # tagE = int(words[0])
                    # dimE = int(words[1])
                    typ = int(words[2])  # type of the element (entity)
                    numElems = int(words[3])  # number of elements of that type
                    for n in range(numElems):
                        l += 1
                        words = lines[l].split()
                        lbl = int(words[0])
                        iNodes = []
                        if (d == 1 and typ == 1):  # 2-node edge
                            iNodes.append(iTab[int(words[1])])
                            iNodes.append(iTab[int(words[2])])
                            mesh.addElement(lbl, typ, iNodes)
                        elif (d == 2 and typ == 2):  # 3-node triangle
                            iNodes.append(iTab[int(words[1])])
                            iNodes.append(iTab[int(words[2])])
                            iNodes.append(iTab[int(words[3])])
                            mesh.addElement(lbl, typ, iNodes)
                        elif (d == 2 and typ == 3):  # 4-node quadrangle
                            iNodes.append(iTab[int(words[1])])
                            iNodes.append(iTab[int(words[2])])
                            iNodes.append(iTab[int(words[3])])
                            iNodes.append(iTab[int(words[4])])
                            mesh.addElement(lbl, typ, iNodes)
            else:
                l += 1  # ignore line
        assert numNode + 1 == nNodes
    else:
        raise ValueError(
            'Format of the meshfile is incorrect. Verify the header $MeshFormat')
    return mesh

#
# Write mesh data in Gmsh v.2 format
# INPUT: - out : output stream
#        - mesh : MeshData object
#


def gmshOutput_mesh(out, mesh):
    """Write mesh in Gmsh v.2 format"""
    dim = mesh.getDimension()
    # HEADER
    out.write("$MeshFormat\n")
    out.write("%d.%d %d %d\n" % (2, 2, 0, 8))
    out.write("$EndMeshFormat\n")
    # NODES
    out.write("$Nodes\n")
    out.write("%d\n" % mesh.nNodes())
    for n in range(mesh.nNodes()):
        out.write("%d" % (n + 1))
        for i in range(dim):
            out.write(" %13.6e" % mesh.getNode(n).getX(i))
        for i in range(dim, 3):
            out.write(" %f" % 0.0)
        out.write("\n")
    out.write("$EndNodes\n")
    # ELEMENTS
    out.write("$Elements\n")
    out.write("%d\n" % mesh.nElements())
    for n in range(mesh.nElements()):
        out.write("%d %d %d %d %d" %
                  (n + 1, mesh.getElement(n).getType(), 2, 1, 1))
        for i in range(mesh.getElement(n).nNodes()):
            out.write(" %d" % (mesh.getElement(n).getNode(i) + 1))
        out.write("\n")
    out.write("$EndElements\n")
    out.flush()

#
# Write nodal field defined on a mesh, in Gmsh V.2 format
# INPUT: - out : output stream
#        - label : string defining the field name
#        - V[nNodes][iDim] : array-like (e.g. NumPy) containing nodal values
#        - it : (time) step number
#        - t : load factor / time
# NOTE: this function must be called after writing the mesh data
#       no attempt is made to check consistency between mesh and field data
#


def gmshOutput_nodal(out, label, V, it, t):
    """Write nodal field in Gmsh format"""
    out.write("$NodeData\n")
    out.write("%d\n" % 1)
    out.write("\"%s\"\n" % label)
    out.write("%d\n" % 1)
    out.write("%f\n" % t)
    out.write("%d\n" % 3)
    out.write("%d\n" % it)
    sz0 = np.shape(V)[1]
    if (sz0 == 2):
        sz1 = 3
    elif(sz0 > 3):
        sz1 = 9
    else:
        sz1 = 1
    out.write("%d\n" % sz1)
    out.write("%d\n" % np.shape(V)[0])
    for n in range(np.shape(V)[0]):
        out.write("%d" % (n + 1))
        for i in range(sz0):
            out.write(" %13.6e" % V[n][i])
        for i in range(sz0, sz1):
            out.write(" %13.6e" % 0.0)
        out.write("\n")
    out.write("$EndNodeData\n")
    out.flush()


#
# Write element field defined on a mesh, in Gmsh V.2 format
# INPUT: - out : output stream
#        - label : string defining the field name
#        - V[nElems][iDim] : array-like (e.g. NumPy) containing element values
#        - it : (time) step number
#        - t : load factor / time
# NOTE: this function must be called after writing the mesh data
#       no attempt is made to check consistency between mesh and field data
#
def gmshOutput_element(out, label, V, it, t):
    """Write element field in Gmsh format"""
    out.write("$ElementData\n")
    out.write("%d\n" % 1)
    out.write("\"%s\"\n" % label)
    out.write("%d\n" % 1)
    out.write("%f\n" % t)
    out.write("%d\n" % 3)
    out.write("%d\n" % it)
    sz0 = np.shape(V)[1]
    if (sz0 == 2 or sz0 == 3):
        sz1 = 3
    elif(sz0 > 3):
        sz1 = 9
    else:
        sz1 = 1
    out.write("%d\n" % sz1)
    out.write("%d\n" % np.shape(V)[0])
    for n in range(np.shape(V)[0]):
        out.write("%d" % (n + 1))
        if (sz1 == 3):
            for i in range(sz0):
                out.write(" %13.6e" % V[n][i])
            for i in range(sz0, sz1):
                out.write(" %13.6e" % 0.0)
        elif (sz1 == 9):
            if (sz0 == 4):
                out.write(" %13.6e %13.6e %13.6e %13.6e %13.6e" %
                          (V[n][0], V[n][1], 0.0, V[n][2], V[n][3]))
                for i in range(4):
                    out.write(" %13.6e" % 0.0)
            else:
                for i in range(sz1):
                    out.write(" %13.6e" % V[n][i])
        out.write("\n")
    out.write("$EndElementData\n")
    out.flush()
