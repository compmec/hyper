# -*- coding: utf-8 -*-
"""
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
    import elasticity
except ModuleNotFoundError as e:
    error = str(e)
    notimportedfilename = error.split("'")[1]
    print("test_residual: Import file not found: " +
          str(notimportedfilename) + ".py")
    sys.exit()

# try:
#     TU.include_path("msh")
# except ModuleNotFoundError:
#     print("Could not include the folder with the mesh examples")


class TUResidual(TU):

    FUNCTIONS = ["new_fileandmaterial", "create_model", "calculate_vectorfunction",
                 "compute_gradstress", "compute_residual", "write_fields", "compare_results"]

    FILENAMES = ["triangle-tri56.msh",
                 "triangle-quad44.msh", "triangle-quad44.msh"]
    MATERIALS = [elasticity.StVenantKirchhoffElasticity,
                 elasticity.NeoHookeanElasticity]

    def __init__(self):
        super(TUResidual, self).__init__()

    def __setUp(self):
        number_files = len(TUResidual.FILENAMES)
        for i in range(number_files):
            TUResidual.FILENAMES[i] = TU.FOLDER_MSH + TUResidual.FILENAMES[i]
        number_materials = len(TUResidual.MATERIALS)

        TUResidual.FUNCTIONS *= number_materials * number_files
        self.set_DEN()
        self.n_mat = 0
        self.n_fil = 0

    def new_fileandmaterial(self):
        if self.n_mat == 0:
            E = 210.e9
            nu = 0.3
        else:
            E = 10.e6
            nu = 0.45
        self.material = TUResidual.MATERIALS[self.n_mat](E, nu)
        self.filename = TUResidual.FILENAMES[self.n_fil]

        # We change the number to make the next test
        self.n_fil += 1
        if self.n_fil >= len(TUResidual.FILENAMES):
            self.n_fil = 0
            self.n_mat += 1

    def create_model(self):
        meshfile = open(self.filename, 'r')
        self.mesh = gmsh.gmshInput_mesh(meshfile)
        meshfile.close()

        # nNodes = self.mesh.getNNodes()
        # nEleme = self.mesh.getNElements()
        # print("Read mesh with %d nodes and %d elements" % (nNodes, nEleme))

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

    def compute_gradstress(self):
        """
        get gradient and stress tensor fields
        """
        # print("self.U = " + str(self.U))
        self.FE_model.computeStress(self.U)

        self.F = []
        self.P = []
        nElements = self.FE_model.getNElements()
        for n in range(nElements):
            new_P = self.FE_model.getStressPK1(n)
            new_F = self.FE_model.getDeformationGradient(n)
            self.F.append(new_F.flatten())
            self.P.append(new_P.flatten())

    def compute_residual(self):
        """
        compute residual array
        """
        # R = np.zeros((mesh.nNodes(), 2))
        # FE_model.computeResidual(U, R)
        self.R = self.FE_model.computeResidual(self.U)

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

    def compare_results(self):
        """
        Not completed, but it should make the comparation between the calculated results (put inside output file) and the reference results (inside reference file).
        For the instant, we don't have a way to comparate the mesh.
        """
        output_filename = self.filename.replace(".msh", "") + '-val.msh'
        reference_filename = self.filename.replace(".msh", "") + '-ref.msh'

        reference_file = open(reference_filename, 'r')
        self.reference_mesh = gmsh.gmshInput_mesh(reference_file)
        reference_file.close()

        output_file = open(output_filename, 'r')
        self.output_mesh = gmsh.gmshInput_mesh(output_file)
        output_file.close()


if __name__ == "__main__":
    test = TUResidual()
    result = test.run()
