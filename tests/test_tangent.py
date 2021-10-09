# -*- coding: utf-8 -*-
"""
# Script to test implementation of the internal tangent matrix in fem.py module

"""

import numpy as np
from hyper import gmsh2 as gmsh
from hyper import elasticity
from hyper import fem


def vector_fct(x):
    fx = x[0] * x[1]
    fy = x[0] * x[1]
    return [fx, fy]


def applyFunctionToMeshNodes(mesh, function):
    nNodes = mesh.getNNodes()
    dim = mesh.getDim()
    U = np.zeros((nNodes, dim))
    for i in range(nNodes):
        Nodei = mesh.getNode(i)
        Xi = Nodei.getX()
        U[i] = vector_fct(Xi)
    return U


def computeKnum(FEmodel, U, DU=1e-5):
    nNodes, dim = U.shape
    Knum = np.zeros((nNodes, dim, nNodes, dim))
    Ud = np.copy(U)

    for i in range(nNodes):
        for j in range(dim):
            Ud[i, j] += DU / 2.
            # R(U_nj + DU/2), residual "plus"
            Rp = FEmodel.computeResidual(Ud)

            Ud[i, j] -= DU
            # R(U_nj - DU/2), residual "minus"
            Rm = FEmodel.computeResidual(Ud)

            # come back to before perturbation for next iterations
            Ud[i, j] += DU / 2
            Knum[i, j] = (Rp - Rm) / DU
    return Knum


def test_Quad11StVenantKirchhoff():
    inputfilename = "../msh/triangle-quad11.msh"
    E = 210.e9
    nu = 0.3
    material = elasticity.StVenantKirchhoffElasticity(E, nu)

    with open(inputfilename, 'r') as meshfile:
        mesh = gmsh.gmshInput_mesh(meshfile)

    FEmodel = fem.FEModel(mesh, material)
    U = applyFunctionToMeshNodes(mesh, vector_fct)
    K = FEmodel.computeTangent(U)
    Knum = computeKnum(FEmodel, U)

    np.testing.assert_allclose(K, Knum)


def test_Quad44StVenantKirchhoff():
    inputfilename = "../msh/triangle-quad44.msh"
    E = 210.e9
    nu = 0.3
    material = elasticity.StVenantKirchhoffElasticity(E, nu)

    with open(inputfilename, 'r') as meshfile:
        mesh = gmsh.gmshInput_mesh(meshfile)

    FEmodel = fem.FEModel(mesh, material)
    U = applyFunctionToMeshNodes(mesh, vector_fct)
    K = FEmodel.computeTangent(U)
    Knum = computeKnum(FEmodel, U)

    np.testing.assert_allclose(K, Knum)


def test_Quad11NeoHookean():
    inputfilename = "../msh/triangle-quad11.msh"
    E = 210.e9
    nu = 0.3
    material = elasticity.StVenantKirchhoffElasticity(E, nu)

    with open(inputfilename, 'r') as meshfile:
        mesh = gmsh.gmshInput_mesh(meshfile)

    FEmodel = fem.FEModel(mesh, material)
    U = applyFunctionToMeshNodes(mesh, vector_fct)
    K = FEmodel.computeTangent(U)
    Knum = computeKnum(FEmodel, U)

    np.testing.assert_allclose(K, Knum)


def test_Quad44NeoHookean():
    inputfilename = "../msh/triangle-quad44.msh"
    E = 210.e9
    nu = 0.3
    material = elasticity.StVenantKirchhoffElasticity(E, nu)

    with open(inputfilename, 'r') as meshfile:
        mesh = gmsh.gmshInput_mesh(meshfile)

    FEmodel = fem.FEModel(mesh, material)
    U = applyFunctionToMeshNodes(mesh, vector_fct)
    K = FEmodel.computeTangent(U)
    Knum = computeKnum(FEmodel, U)

    np.testing.assert_allclose(K, Knum)
