# -*- coding: utf-8 -*-
#
# Script to test implementation of the discrete residual
#

"""
MNMNL Homework
Module 5: discretized residual

This script tests your implementation of the discretized residual in fem.py module

You will test two types of constitutive models:
    - StVenantKirchhoff compressible material
    - NeoHookean compressible material

You will test two types of meshes:
    - triangle-tri56.msh with 56 triangular elements
    - triangle-quad44.msh with 44 quadrangular elements

Compare your residual field with the solution given in the mesh files:
    - triangle-tri56-ref.msh with 56 triangular elements
    - triangle-quad44-ref.msh with 44 quadrangular elements
compare_residual.geo script, that must be opened with Gmsh.

WARNING: the reference solutions are only computed for the NeoHookean model !!!
"""
#
# --- namespace
#
# import numpy to create arrays (access numpy methods with np.method())
# import elasticity.py (from previous homework)

# import completed fem.py (to complete)


try:
    import numpy as np
except ModuleNotFoundError:
    print("The numpy library was not found")


try:
    from test_unit import TU
except ModuleNotFoundError:
    print("Main file to test was not found: test_unit.py")

try:
    TU.include_path("src")
    import gmsh2 as gmsh
    import elasticity
    import fem
except ModuleNotFoundError:
    print("Import file not found: gmsh2.py")


class TUResidual(TU):

    FUNCTIONS = ["calculate_vectorfunction", "compute_residual",
                 "compute_gradstress", "write_fields"]

    def __init__(self, inputfilename, material):
        self.filename = inputfilename
        self.material = material
        super(TUResidual, self).__init__()

    def __setUp(self):
        meshfile = open(self.filename, 'r')
        self.mesh = gmsh.gmshInput_mesh(meshfile)
        meshfile.close()

        nNodes = self.mesh.getNNodes()
        nEleme = self.mesh.getNElements()
        print("Read mesh with %d nodes and %d elements" % (nNodes, nEleme))

        self.FE_model = fem.FEModel(self.mesh, self.material)

    def calculate_vectorfunction(self):
        def fct(x):
            fx = x[0] * x[1]
            fy = x[0] * x[1]
            return [fx, fy]

        #
        # -- create nodal array
        #
        nNodes = self.mesh.getNNodes()
        dim = self.mesh.getDim()
        self.U = np.zeros((nNodes, dim))  # initialise array of zeros,
        # of shape nNodes x 2: [U1_x,U1_y,...,UnNodes_x,UnNodes_y]
        for i in range(nNodes):
            Nodei = self.mesh.getNode(i)
            Xi = Nodei.getX()
            self.U[i] = fct(Xi)

    def compute_residual(self):
        #
        # -- compute residual array
        #
        # R = np.zeros((mesh.nNodes(), 2))
        # FE_model.computeResidual(U, R)
        self.R = self.FE_model.computeResidual(self.U)

    def compute_gradstress(self):
        #
        # -- get gradient and stress tensor fields
        #
        self.FE_model.computeStress(self.U)

        self.F = []
        self.P = []
        nElements = self.FE_model.getNElements()
        for n in range(nElements):
            new_F = self.FE_model.getDeformationGradient(n)
            new_P = self.FE_model.getStressPK1(n)
            self.F.append(new_F.flatten())
            self.P.append(new_P.flatten())

    def write_fields(self):
        outputfnam = self.filename.replace(".msh", "") + '-val.msh'
        outfile = open(outputfnam, 'w')
        gmsh.gmshOutput_mesh(outfile, self.mesh)
        gmsh.gmshOutput_nodal(outfile, "displacement", self.U, 0, 0.0)
        gmsh.gmshOutput_nodal(outfile, "residual", self.R, 0, 0.0)
        gmsh.gmshOutput_element(
            outfile, "deformation gradient", self.F, 0, 0.0)
        gmsh.gmshOutput_element(outfile, "engineering stress", self.P, 0, 0.0)
        outfile.close()


if __name__ == "__main__":
    msh_folder = "../msh/"
    files = ["triangle-tri56.msh", "triangle-quad44.msh"]
    for i in range(2):
        if i == 0:
            E = 210.e9
            nu = 0.3
            mat_model = elasticity.StVenantKirchhoffElasticity(E, nu)
        elif i == 1:
            E = 10.e6
            nu = 0.45
            mat_model = elasticity.NeoHookeanElasticity(E, nu)
        for file in files:
            inputfilename = msh_folder + file
            test = TUResidual(inputfilename, mat_model)
            test.run()
