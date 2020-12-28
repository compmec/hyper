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
except ModuleNotFoundError:
    print("Import file not found: gmsh2.py")

try:
    TU.include_path("msh")
except ModuleNotFoundError:
    print("Could not include the folder with the mesh examples")


def fct(x):
    f = x[0]**2 + x[1]**2
    return f


def vector_fct(x):
    fx = 2 * x[0]
    fy = -0.2 * x[1]
    return np.array((fx, fy))


class TUMesh(TU):

    FUNCTIONS = ["calculate_scalarfunction", "calculate_vectorfield",
                 "write_mesh", "write_scalarfunction", "write_vectorfield"]

    def __init__(self, inputfilename):
        super(TUMesh, self).__init__()
        self.filename = inputfilename

    def __setUp(self):
        meshfile = open(self.filename, 'r')
        self.mesh = gmsh.gmshInput_mesh(meshfile)
        meshfile.close()
        print("Read mesh with %d nodes and %d elements" %
              (self.mesh.nNodes(), self.mesh.nElements()))

    def calculate_scalarfunction(self):
        # initialise array of zeros,
        nNodes = self.mesh.nNodes()
        self.V = np.zeros((nNodes, 1))
        for i in range(self.mesh.nNodes()):
            Node_i = self.mesh.getNode(i)
            X_i = Node_i.getX()
            self.V[i] = fct(X_i)

    def calculate_vectorfield(self):
        # --- create nodal array
        nNodes = self.mesh.nNodes()
        dim = self.mesh.getDim()
        self.U = np.zeros((nNodes, dim))  # initialise array of zeros,
        # We have that U is a array of Nx2.
        # That is, U = [[U0x, U0y], [U1x, U1y], ..., [UNx, UNy]]
        # If we want U to be just a vector
        # U = U.reshape(U.size)

        # print("Calculating vector field for %d nodes" % (self.mesh.nNodes()))
        for i in range(nNodes):
            Node_i = self.mesh.getNode(i)
            X_i = Node_i.getX()
            self.U[i] = vector_fct(X_i)

    def write_mesh(self):
        outputfnam = self.filename.replace(".msh", "") + '-val.msh'
        self.outfile = open(outputfnam, 'w')
        # print("Writing mesh into " + str(outputfnam))
        gmsh.gmshOutput_mesh(self.outfile, self.mesh)

    def write_scalarfunction(self):
        # print("Writing scalar field into " + str(outputfnam))
        gmsh.gmshOutput_nodal(self.outfile, "scalar field", self.V,
                              it=0, t=0.0)  # write nodal field 'V'
        # named 'scalar field' to self.mesh 'outfile' at iteration 'it' and time 't'

    def write_vectorfield(self):
        # print("Writing vector field in " + str(outputfnam))
        gmsh.gmshOutput_nodal(self.outfile, "VectorField", self.U, it=0, t=0)
        self.outfile.close()


if __name__ == "__main__":
    msh_folder = "../msh/"
    files = ["meshfile-tri8.msh", "meshfile-quad8.msh"]

    N = len(files)  # Number of tests
    for i in range(N):
        print("##########################")
        print("#   Test File Number " + str(i + 1) + "   #")
        print("##########################")
        filename = msh_folder + files[i]
        test = TUMesh(filename)
        test.run()
