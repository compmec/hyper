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

from hyper import fem
from hyper import gmsh2 as gmsh


def vector_fct(x):
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


def test_Tri14():
    inputfilename = "tests/msh/triangle-tri14.msh"
    with open(inputfilename, "r") as meshfile:
        mesh = gmsh.gmshInput_mesh(meshfile)
    nNodes = mesh.nNodes()
    Nu = mesh.getDim()  # Number of unknows by for each point
    U = np.zeros((nNodes, Nu))  # initializes array of zeros,
    for i in range(mesh.nNodes()):
        Nodei = mesh.getNode(i)
        Xi = Nodei.getX()
        U[i] = vector_fct(Xi)

    # Create FE model
    FE_model = fem.FEModel(mesh)

    # Compute deformation gradient and langrangian strain tensor fields
    FE_model.computeStrain(U)  # compute gradient and strain
    F = []  # list of deformation gradient by element
    E = []  # list of lagrangian strain by element
    nElements = FE_model.nElements()
    for n in range(nElements):
        F.append(FE_model.getDeformationGradient(n).flatten())
        E.append(FE_model.getStrainGreenLagrange(n).flatten())

    outputfilename = inputfilename.replace(".msh", "-val.msh")
    with open(outputfilename, "w") as outputfile:
        gmsh.gmshOutput_mesh(outputfile, mesh)

        fieldname = "U: displacement field"
        gmsh.gmshOutput_nodal(outputfile, fieldname, U, it=0, t=0)

        fieldname = "F: Deformation Gradient"
        gmsh.gmshOutput_element(outputfile, fieldname, F, it=0, t=0)

        fieldname = "E: Strain Green Lagrange"
        gmsh.gmshOutput_element(outputfile, fieldname, E, it=0, t=0)


def test_Tri56():
    inputfilename = "tests/msh/triangle-tri56.msh"
    with open(inputfilename, "r") as meshfile:
        mesh = gmsh.gmshInput_mesh(meshfile)
    nNodes = mesh.nNodes()
    Nu = mesh.getDim()  # Number of unknows by for each point
    U = np.zeros((nNodes, Nu))  # initializes array of zeros,
    for i in range(mesh.nNodes()):
        Nodei = mesh.getNode(i)
        Xi = Nodei.getX()
        U[i] = vector_fct(Xi)

    # Create FE model
    FE_model = fem.FEModel(mesh)

    # Compute deformation gradient and langrangian strain tensor fields
    FE_model.computeStrain(U)  # compute gradient and strain
    F = []  # list of deformation gradient by element
    E = []  # list of lagrangian strain by element
    nElements = FE_model.nElements()
    for n in range(nElements):
        F.append(FE_model.getDeformationGradient(n).flatten())
        E.append(FE_model.getStrainGreenLagrange(n).flatten())

    outputfilename = inputfilename.replace(".msh", "-val.msh")
    with open(outputfilename, "w") as outputfile:
        gmsh.gmshOutput_mesh(outputfile, mesh)

        fieldname = "U: displacement field"
        gmsh.gmshOutput_nodal(outputfile, fieldname, U, it=0, t=0)

        fieldname = "F: Deformation Gradient"
        gmsh.gmshOutput_element(outputfile, fieldname, F, it=0, t=0)

        fieldname = "E: Strain Green Lagrange"
        gmsh.gmshOutput_element(outputfile, fieldname, E, it=0, t=0)


def test_Quad11():
    inputfilename = "tests/msh/triangle-quad11.msh"
    with open(inputfilename, "r") as meshfile:
        mesh = gmsh.gmshInput_mesh(meshfile)
    nNodes = mesh.nNodes()
    Nu = mesh.getDim()  # Number of unknows by for each point
    U = np.zeros((nNodes, Nu))  # initializes array of zeros,
    for i in range(mesh.nNodes()):
        Nodei = mesh.getNode(i)
        Xi = Nodei.getX()
        U[i] = vector_fct(Xi)

    # Create FE model
    FE_model = fem.FEModel(mesh)

    # Compute deformation gradient and langrangian strain tensor fields
    FE_model.computeStrain(U)  # compute gradient and strain
    F = []  # list of deformation gradient by element
    E = []  # list of lagrangian strain by element
    nElements = FE_model.nElements()
    for n in range(nElements):
        F.append(FE_model.getDeformationGradient(n).flatten())
        E.append(FE_model.getStrainGreenLagrange(n).flatten())

    outputfilename = inputfilename.replace(".msh", "-val.msh")
    with open(outputfilename, "w") as outputfile:
        gmsh.gmshOutput_mesh(outputfile, mesh)

        fieldname = "U: displacement field"
        gmsh.gmshOutput_nodal(outputfile, fieldname, U, it=0, t=0)

        fieldname = "F: Deformation Gradient"
        gmsh.gmshOutput_element(outputfile, fieldname, F, it=0, t=0)

        fieldname = "E: Strain Green Lagrange"
        gmsh.gmshOutput_element(outputfile, fieldname, E, it=0, t=0)


def test_Quad44():
    inputfilename = "tests/msh/triangle-quad44.msh"
    with open(inputfilename, "r") as meshfile:
        mesh = gmsh.gmshInput_mesh(meshfile)
    nNodes = mesh.nNodes()
    Nu = mesh.getDim()  # Number of unknows by for each point
    U = np.zeros((nNodes, Nu))  # initializes array of zeros,
    for i in range(mesh.nNodes()):
        Nodei = mesh.getNode(i)
        Xi = Nodei.getX()
        U[i] = vector_fct(Xi)

    # Create FE model
    FE_model = fem.FEModel(mesh)

    # Compute deformation gradient and langrangian strain tensor fields
    FE_model.computeStrain(U)  # compute gradient and strain
    F = []  # list of deformation gradient by element
    E = []  # list of lagrangian strain by element
    nElements = FE_model.nElements()
    for n in range(nElements):
        F.append(FE_model.getDeformationGradient(n).flatten())
        E.append(FE_model.getStrainGreenLagrange(n).flatten())

    outputfilename = inputfilename.replace(".msh", "-val.msh")
    with open(outputfilename, "w") as outputfile:
        gmsh.gmshOutput_mesh(outputfile, mesh)

        fieldname = "U: displacement field"
        gmsh.gmshOutput_nodal(outputfile, fieldname, U, it=0, t=0)

        fieldname = "F: Deformation Gradient"
        gmsh.gmshOutput_element(outputfile, fieldname, F, it=0, t=0)

        fieldname = "E: Strain Green Lagrange"
        gmsh.gmshOutput_element(outputfile, fieldname, E, it=0, t=0)
