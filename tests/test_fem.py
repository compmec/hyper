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


try:
    import sys
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
    import fem
except ModuleNotFoundError as e:
    error = str(e)
    notimportedfilename = error.split("'")[1]
    print("Import file not found: " + str(notimportedfilename) + ".py")
    sys.exit()


class TUFem(TU):

    FUNCTIONS = ["create_nodalarray", "compute_deformation",
                 "write_scalarfields", "write_tensorfields"]

    def __init__(self, inputfilename):
        self.filename = inputfilename
        super(TUFem, self).__init__()

    def __setUp(self):
        meshfile = open(self.filename, 'r')
        self.mesh = gmsh.gmshInput_mesh(meshfile)
        meshfile.close()
        print("Read mesh with %d nodes and %d elements" %
              (self.mesh.nNodes(), self.mesh.nElements()))

    def create_nodalarray(self):
        #
        # --- create nodal array
        #
        nNodes = self.mesh.nNodes()
        Nu = self.mesh.getDim()  # Number of unknows by for each point
        self.U = np.zeros((nNodes, Nu))  # initializes array of zeros,
        for i in range(self.mesh.nNodes()):
            Nodei = self.mesh.getNode(i)
            Xi = Nodei.getX()
            self.U[i] = fct(Xi)

    def compute_deformation(self):
        # --- create FE model
        FE_model = fem.FEModel(self.mesh)

        # --- compute deformation gradient and langrangian strain tensor fields
        FE_model.computeStrain(self.U)  # compute gradient and strain
        self.F = []  # list of deformation gradient by element
        self.E = []  # list of lagrangian strain by element
        for n in range(FE_model.nElements()):
            # flatten() transforms a 2D array into a 1D array, cf. numpy documentation
            self.F.append(FE_model.getDeformationGradient(n).flatten())
            # flatten() transforms a 2D array into a 1D array, cf. numpy documentation
            self.E.append(FE_model.getStrainGreenLagrange(n).flatten())

    def write_scalarfields(self):
        # --- write displacement field
        # write mesh to the file (only do this once!)
        outputfilename = self.filename.replace(".msh", "") + '-val.msh'
        self.outfile = open(outputfilename, 'w')  # open new meshfile
        gmsh.gmshOutput_mesh(self.outfile, self.mesh)
        gmsh.gmshOutput_nodal(
            self.outfile, "U: displacement field", self.U, 0, 0.0)

    def write_tensorfields(self):
        # print("Writing F in " + str(outputfilename))
        gmsh.gmshOutput_element(self.outfile, "F", self.F, it=0, t=0)
        # print("Writing E in " + str(outputfilename))
        gmsh.gmshOutput_element(self.outfile, "E", self.E, it=0, t=0)
        self.outfile.close()


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

    files = ["triangle-tri14", "triangle-tri56",
             "triangle-quad11", "triangle-quad44"]
    for basefile in files:
        inputfilename = msh_folder + basefile + ".msh"
        test = TUFem(inputfilename)
        test.run()

    #
    # --- open Gmsh output file
    #

    #
    # --- write tensor fields to gmsh output file
    #
