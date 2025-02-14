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

import hyper.gmsh2 as gmsh


def scalar_fct(x):
    f = x[0] ** 2 + x[1] ** 2
    return f


def vector_fct(x):
    fx = 2 * x[0]
    fy = -0.2 * x[1]
    return np.array((fx, fy))


def test_ReadTri8():
    # Reading the mesh
    inputfilename = "tests/msh/meshfile-tri8.msh"
    with open(inputfilename, "r") as meshfile:
        mesh = gmsh.gmshInput_mesh(meshfile)

    # Testing the number of points and dimension
    assert mesh.nNodes() == 8
    assert mesh.getDim() == 2
    assert mesh.getNElements() == 8


def test_ReadTri14():
    # Reading the mesh
    inputfilename = "tests/msh/triangle-tri14.msh"
    with open(inputfilename, "r") as meshfile:
        mesh = gmsh.gmshInput_mesh(meshfile)

    # Testing the number of points and dimension
    assert mesh.nNodes() == 13
    assert mesh.getDim() == 2
    assert mesh.getNElements() == 14


def test_ReadTri56():
    # Reading the mesh
    inputfilename = "tests/msh/triangle-tri56.msh"
    with open(inputfilename, "r") as meshfile:
        mesh = gmsh.gmshInput_mesh(meshfile)

    assert mesh.nNodes() == 39
    assert mesh.getDim() == 2
    assert mesh.getNElements() == 56


def test_ReadQuad8():
    # Reading the mesh
    inputfilename = "tests/msh/meshfile-quad8.msh"
    with open(inputfilename, "r") as meshfile:
        mesh = gmsh.gmshInput_mesh(meshfile)

    # Testing the number of points and dimension
    assert mesh.nNodes() == 15
    assert mesh.getDim() == 2
    assert mesh.getNElements() == 8


def test_ReadQuad11():
    # Reading the mesh
    inputfilename = "tests/msh/triangle-quad11.msh"
    with open(inputfilename, "r") as meshfile:
        mesh = gmsh.gmshInput_mesh(meshfile)

    # Testing the number of points and dimension
    assert mesh.nNodes() == 16
    assert mesh.getDim() == 2
    assert mesh.getNElements() == 11


def test_ReadQuad44():
    # Reading the mesh
    inputfilename = "tests/msh/triangle-quad44.msh"
    with open(inputfilename, "r") as meshfile:
        mesh = gmsh.gmshInput_mesh(meshfile)

    # Testing the number of points and dimension
    assert mesh.nNodes() == 49
    assert mesh.getDim() == 2
    assert mesh.getNElements() == 44


############################################################
#    TESTING THE POSITION OF NODES
#################################################


def test_PositionNodesTri8():
    # Reading the mesh
    inputfilename = "tests/msh/meshfile-tri8.msh"
    with open(inputfilename, "r") as meshfile:
        mesh = gmsh.gmshInput_mesh(meshfile)

    # Testing the position of the node
    nNodes = mesh.nNodes()
    dim = mesh.getDim()
    X = np.zeros((nNodes, dim))
    for i in range(nNodes):
        Node_i = mesh.getNode(i)
        X[i] = Node_i.getX()
    Xgood = [
        [0, 1, 1, 0, 0.5, 0.5, 0.25, 0.75],
        [0, 0, 0.3, 0.3, 0, 0.3, 0.15, 0.15],
    ]
    Xgood = np.array(Xgood).T
    np.testing.assert_almost_equal(X, Xgood)


def test_PositionNodesQuad8():
    inputfilename = "tests/msh/meshfile-quad8.msh"
    with open(inputfilename, "r") as meshfile:
        mesh = gmsh.gmshInput_mesh(meshfile)

    # Testing the position of the node
    nNodes = mesh.nNodes()
    dim = mesh.getDim()
    X = np.zeros((nNodes, dim))
    for i in range(nNodes):
        Node_i = mesh.getNode(i)
        X[i] = Node_i.getX()
    Xgood = [
        [0, 20, 20, 0, 5, 10, 15, 20, 15, 10, 5, 0, 10, 5, 15],
        [0, 0, 6, 6, 0, 0, 0, 3, 6, 6, 6, 3, 3, 3, 3],
    ]
    Xgood = np.array(Xgood).T / 20
    np.testing.assert_almost_equal(X, Xgood)


# ################################################
# #  TESTING FIELDS
# ################################################


def test_FieldsTri8():
    inputfilename = "tests/msh/meshfile-tri8.msh"
    with open(inputfilename, "r") as meshfile:
        mesh = gmsh.gmshInput_mesh(meshfile)
    nNodes = mesh.nNodes()
    dim = mesh.getDim()

    # Testing the scalar field
    V = np.zeros((nNodes, 1))
    for i in range(nNodes):
        Node_i = mesh.getNode(i)
        X_i = Node_i.getX()
        V[i] = scalar_fct(X_i)

    # Testing the vectorial field
    U = np.zeros((nNodes, dim))
    for i in range(nNodes):
        Node_i = mesh.getNode(i)
        X_i = Node_i.getX()
        U[i] = vector_fct(X_i)

    # Writting the results in a mesh
    outputfilename = inputfilename.replace(".msh", "-val.msh")
    with open(outputfilename, "w") as outputfile:
        # Write the base mesh
        gmsh.gmshOutput_mesh(outputfile, mesh)

        # Write scalar field
        gmsh.gmshOutput_nodal(
            outputfile, "scalar field", V, it=0, t=0.0
        )  # write nodal field 'V'

        # Write vectorial field
        gmsh.gmshOutput_nodal(outputfile, "VectorField", U, it=0, t=0)


def test_FieldsQuad8():
    inputfilename = "tests/msh/meshfile-quad8.msh"
    with open(inputfilename, "r") as meshfile:
        mesh = gmsh.gmshInput_mesh(meshfile)
    nNodes = mesh.nNodes()
    dim = mesh.getDim()

    # Testing the scalar field
    V = np.zeros((nNodes, 1))
    for i in range(nNodes):
        Node_i = mesh.getNode(i)
        X_i = Node_i.getX()
        V[i] = scalar_fct(X_i)

    # Testing the vectorial field
    U = np.zeros((nNodes, dim))
    for i in range(nNodes):
        Node_i = mesh.getNode(i)
        X_i = Node_i.getX()
        U[i] = vector_fct(X_i)

    # Writting the results in a mesh
    outputfilename = inputfilename.replace(".msh", "-val.msh")
    with open(outputfilename, "w") as outputfile:
        # Write the base mesh
        gmsh.gmshOutput_mesh(outputfile, mesh)

        # Write scalar field
        gmsh.gmshOutput_nodal(
            outputfile, "scalar field", V, it=0, t=0.0
        )  # write nodal field 'V'

        # Write vectorial field
        gmsh.gmshOutput_nodal(outputfile, "VectorField", U, it=0, t=0)
