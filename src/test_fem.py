# -*- coding: utf-8 -*-
#
# Script to test kinematics of fem module
#
"""
MNMNL Homework
Module 3: kinematic tensor field on a mesh
complete code below where indicated by
# ============
# TODO:
# ...
# ============

Test your implementation of the fem.py module which uses tensor.py for the problem
given in the exercise.
You will test two types of meshes (with two refinements):
    - triangle-tri(14,56).msh with 14 and 56 triangular elements respectively
    - meshfile-quad(11,44).msh with 11 and 44 quadrangular elements respectively
Compare your nodal and element fields with the anaylitical solution you have
found in the exercise.
"""


import numpy as np
import gmsh
import fem


def fct(x):
    """
    Function that translates the position X into the vector U displacement
    Field:
        U(X1,X2) = X1*X2(a*e1 + b*e2)
    # for a=b=1 (Second and Third part of the exercise)
    """
    a = b = 1
    f0 = a * x[0] * x[1]
    f1 = b * x[0] * x[1]
    return [f0, f1]


if __name__ == "__main__":

    #
    # --- read mesh from Gmsh file
    #
    msh_folder = "../msh/"

    basenamefile = "triangle-tri14"
    # basenamefile = "triangle-tri56"
    # basenamefile = "triangle-quad11"
    # basenamefile = "triangle-quad44"

    inputfnam = msh_folder + basenamefile + ".msh"
    meshfile = open(inputfnam, 'r')
    mesh = gmsh.gmshInput_mesh(meshfile)
    meshfile.close()
    print("Read mesh with %d nodes and %d elements" %
          (mesh.nNodes(), mesh.nElements()))

    #
    # --- open Gmsh output file
    #
    outputfilename = msh_folder + basenamefile + '-val.msh'
    outfile = open(outputfilename, 'w')  # open new meshfile

    #
    # --- create nodal array
    #
    nNodes = mesh.nNodes()
    Nu = 2  # Number of unknows by for each point
    U = np.zeros((nNodes, Nu))  # initializes array of zeros,
    for i in range(mesh.nNodes()):
        Xi = mesh.getNode(i).getX()
        U[i] = fct(Xi)

    #
    # --- write displacement field
    #
    # write mesh to the file (only do this once!)
    gmsh.gmshOutput_mesh(outfile, mesh)
    gmsh.gmshOutput_nodal(outfile, "U: displacement field", U, 0, 0.0)

    #
    # --- create FE model
    #
    FE_model = fem.FEModel(mesh)

    #
    # --- compute deformation gradient and langrangian strain tensor fields
    #
    FE_model.computeStrain(U)  # compute gradient and strain
    F = []  # list of deformation gradient by element
    E = []  # list of lagrangian strain by element
    for n in range(FE_model.nElements()):
        # flatten() transforms a 2D array into a 1D array, cf. numpy documentation
        F.append(FE_model.getDeformationGradient(n).flatten())
        # flatten() transforms a 2D array into a 1D array, cf. numpy documentation
        E.append(FE_model.getStrainGreenLagrange(n).flatten())

    #
    # --- write tensor fields to gmsh output file
    #
    print("Writing F in " + str(outputfilename))
    gmsh.gmshOutput_element(outfile, "F", F, it=0, t=0)
    print("Writing E in " + str(outputfilename))
    gmsh.gmshOutput_element(outfile, "E", E, it=0, t=0)
    # =============================================================================

    #
    # --- close Gmsh output mesh file after all input/output operations
    #
    outfile.close()
