# -*- coding: utf-8 -*-
"""
"""


import numpy as np

from hyper import elasticity, fem
from hyper import gmsh2 as gmsh


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


def computeTestValues(FEmodel, U):
    FEmodel.computeStress(U)
    F = []
    P = []
    nElements = FEmodel.getNElements()
    for n in range(nElements):
        new_P = FEmodel.getStressPK1(n)
        new_F = FEmodel.getDeformationGradient(n)
        F.append(new_F.flatten())
        P.append(new_P.flatten())
    R = FEmodel.computeResidual(U)
    return R, F, P


def writeFieldsOnGmshFile(fields, mesh, filename):
    with open(filename, "w") as file:
        gmsh.gmshOutput_mesh(file, mesh)

        for fieldname, field in fields["nodal"]:
            gmsh.gmshOutput_nodal(file, fieldname, field, 0, 0.0)

        for fieldname, field in fields["element"]:
            gmsh.gmshOutput_element(file, fieldname, field, 0, 0.0)


def test_Tri56StVenantKirchhoff():
    inputfilename = "tests/msh/triangle-tri56.msh"
    outputfilename = "tests/msh/triangle-tri56-StVenant-val.msh"
    E = 210.0e9
    nu = 0.3
    material = elasticity.StVenantKirchhoffElasticity(E, nu)

    with open(inputfilename, "r") as meshfile:
        mesh = gmsh.gmshInput_mesh(meshfile)

    FEmodel = fem.FEModel(mesh, material)
    U = applyFunctionToMeshNodes(mesh, vector_fct)

    R, F, P = computeTestValues(FEmodel, U)

    fields = {
        "nodal": [("U: Displacement", U), ("R: Residual", R)],
        "element": [
            ("F: Deformation Gradient", F),
            ("P: Engineering stress", P),
        ],
    }
    writeFieldsOnGmshFile(fields, mesh, outputfilename)


def test_Quad44StVenantKirchhoff():
    inputfilename = "tests/msh/triangle-quad44.msh"
    outputfilename = "tests/msh/triangle-quad44-StVenant-val.msh"
    E = 210.0e9
    nu = 0.3
    material = elasticity.StVenantKirchhoffElasticity(E, nu)

    with open(inputfilename, "r") as meshfile:
        mesh = gmsh.gmshInput_mesh(meshfile)

    FEmodel = fem.FEModel(mesh, material)
    U = applyFunctionToMeshNodes(mesh, vector_fct)

    R, F, P = computeTestValues(FEmodel, U)

    fields = {
        "nodal": [("U: Displacement", U), ("R: Residual", R)],
        "element": [
            ("F: Deformation Gradient", F),
            ("P: Engineering stress", P),
        ],
    }
    writeFieldsOnGmshFile(fields, mesh, outputfilename)


def test_Tri56NeoHookean():
    inputfilename = "tests/msh/triangle-tri56.msh"
    outputfilename = "tests/msh/triangle-tri56-NeoHookean-val.msh"
    E = 10.0e6
    nu = 0.45
    material = elasticity.StVenantKirchhoffElasticity(E, nu)

    with open(inputfilename, "r") as meshfile:
        mesh = gmsh.gmshInput_mesh(meshfile)

    FEmodel = fem.FEModel(mesh, material)
    U = applyFunctionToMeshNodes(mesh, vector_fct)

    R, F, P = computeTestValues(FEmodel, U)

    fields = {
        "nodal": [("U: Displacement", U), ("R: Residual", R)],
        "element": [
            ("F: Deformation Gradient", F),
            ("P: Engineering stress", P),
        ],
    }
    writeFieldsOnGmshFile(fields, mesh, outputfilename)


def test_Quad44NeoHookean():
    inputfilename = "tests/msh/triangle-tri56.msh"
    outputfilename = "tests/msh/triangle-tri56-NeoHookean-val.msh"
    E = 10.0e6
    nu = 0.45
    material = elasticity.StVenantKirchhoffElasticity(E, nu)

    with open(inputfilename, "r") as meshfile:
        mesh = gmsh.gmshInput_mesh(meshfile)

    FEmodel = fem.FEModel(mesh, material)
    U = applyFunctionToMeshNodes(mesh, vector_fct)

    R, F, P = computeTestValues(FEmodel, U)

    fields = {
        "nodal": [("U: Displacement", U), ("R: Residual", R)],
        "element": [
            ("F: Deformation Gradient", F),
            ("P: Engineering stress", P),
        ],
    }
    writeFieldsOnGmshFile(fields, mesh, outputfilename)


# def compare_results(self):
#     """
#     Not completed, but it should make the comparation between the calculated results (put inside output file) and the reference results (inside reference file).
#     For the instant, we don't have a way to comparate the mesh.
#     """
#     output_filename = self.filename.replace(".msh", "") + '-val.msh'
#     reference_filename = self.filename.replace(".msh", "") + '-ref.msh'

#     reference_file = open(reference_filename, 'r')
#     self.reference_mesh = gmsh.gmshInput_mesh(reference_file)
#     reference_file.close()

#     output_file = open(output_filename, 'r')
#     self.output_mesh = gmsh.gmshInput_mesh(output_file)
#     output_file.close()


# Mesh libraries
# open3d, plyfile, pymesh, pygmsh, meshio
# https://newbedev.com/python-plyfile-vs-pymesh
