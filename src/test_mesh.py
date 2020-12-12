# -*- coding: utf-8 -*-
#
# Script to test my_mesh module
#
"""
MNMNL Homework
Module 2: nodal field on a my_mesh complete code below where indicated by
    Carlos Adir ELY MURUSSI LEITE
    Clément FUCHS
# ============
Test your implementation of the my_mesh.py module on two different meshes:
    - meshfile-tri8.msh with 8 triangular elements
    - meshfile-quad8.msh with 8 quadrangular elements
Compare your nodal fields with the solution given in the my_mesh files:
    - meshfile-tri8-ref.msh with 8 triangular elements
    - meshfile-quad8-ref.msh with 8 quadrangular elements
compare_meshfile.geo script, that must be opened with Gmsh.
"""

import numpy as np
import gmsh


def fct(x):
    f = x[0]**2 + x[1]**2
    return f


def vector_fct(x):
    fx = 2 * x[0]
    fy = -0.2 * x[1]
    return np.array((fx, fy))


if __name__ == "__main__":

    geo_folder = "../geo/"
    msh_folder = "../msh/"

    # inputfnam = msh_folder + "meshfile-tri8.msh"
    inputfnam = msh_folder + "meshfile-quad8.msh"
    meshfile = open(inputfnam, 'r')
    my_mesh = gmsh.gmshInput_mesh(meshfile)
    meshfile.close()

    print("Read mesh with %d nodes and %d elements" %
          (my_mesh.nNodes(), my_mesh.nElements()))

    # =========================================================================
    #  Scalar field
    # =========================================================================
    # --- create nodal array
    V = np.zeros((my_mesh.nNodes(), 1))  # initialise array of zeros,
    # Here we create a numpy matrix of Nx1 instead of a numpy vector of N
    # elements. The reason is that in the memory they always stores as matrix
    # So, it doesn't matter very much.
    print("Calculating scalar field for %d nodes" % (my_mesh.nNodes()))
    for i in range(my_mesh.nNodes()):
        Node_i = my_mesh.getNode(i)
        X_i = Node_i.getX()
        V[i] = fct(X_i)

    # --- write scalar field
    outputfnam = inputfnam.replace(".msh", "") + '-val.msh'
    outfile = open(outputfnam, 'w')
    print("Writing mesh into " + str(outputfnam))
    gmsh.gmshOutput_mesh(outfile, my_mesh)
    print("Writing scalar field into " + str(outputfnam))
    gmsh.gmshOutput_nodal(outfile, "scalar field", V, it=0,
                          t=0.0)  # write nodal field 'V'
    # named 'scalar field' to my_mesh 'outfile' at iteration 'it' and time 't'

    # =========================================================================
    # vector field
    # =========================================================================
    # --- create nodal array
    U = np.zeros((my_mesh.nNodes(), 2))  # initialise array of zeros,
    # We have that U is a array of Nx2.
    # That is, U = [[U0x, U0y], [U1x, U1y], ..., [UNx, UNy]]
    # If we want U to be just a vector
    # U = U.reshape(U.size)

    print("Calculating vector field for %d nodes" % (my_mesh.nNodes()))
    for i in range(my_mesh.nNodes()):
        Node_i = my_mesh.getNode(i)
        X_i = Node_i.getX()
        U[i] = vector_fct(X_i)

    # --- write vector field
    print("Writing vector field in " + str(outputfnam))
    gmsh.gmshOutput_nodal(outfile, "VectorField", U, it=0, t=0)

    # --- close my_mesh file after all input/output operations
    outfile.close()
