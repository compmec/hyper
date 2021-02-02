# -*- coding: utf-8 -*-
#
# Script to test implementation of the internal tangent matrix in fem.py module
#
"""
MNMNL Homework
Module 6: discretized tangent matrix

This script tests your implementation of the discretized tangent matrix in fem.py
module with a numerical differentiation of the discretized residual with respect
to the displacements.

You MUST test both material models:
    - NeoHookean
    - St Venant Kirchhoff

You can test whichever meshfile you want.

WARNING: THERE IS NOTHING TO COMPLETE HERE !

TO DO:
1) You must complete fem.py module as follows:
    - FiniteElement.computeStiffness():
        computes Knod, internal tangent matrix of the element e, shape: (nNodes(e), dim, nNodes(e), dim)
    - FEModel.assemble4():
        can be used to assemble tangent matrix and called in computeInternalStiffness
    - FEModel.computeInternalStiffness():
        computes K, assembled internal tangent matrix of the whole mesh, shape: (nNodes(total),dim,nNodes(total),dim)
2) You must complete tensor.py as follows:
    - MaterialToLagrangian(F,S,M):
        computes dP/dF, the Lagrangian tangent operator from F the deformation gradient tensor,
        S the second Piola-Kirchhoff stress tensor, and M = 2*dS/dC the material tangent operator.

WARNING 2:
in FEModel.computeInternalStiffness(), K is an argument and is modified directly in the function.
Therefore, the first line
K[:,:,:,:] = 0.0
sets all components of K to zero : K is erased at each call of the function and computed again.
It saves memory (no need to create an empty array (hence allocate memory for it) at each call)

WARNING 3:
To validate the code, your computation of the residual must be correct. Be sure to
correct the M5 assignment if necessary !
"""

#
#--- namespace
#
# import numpy to create arrays (access numpy methods with np.method())
import sys
sys.path.append("../src/")
import numpy as np
# import gmsh.py module (given)
import gmsh2 as gmsh
# import elasticity.py (from previous homework)
import elasticity
# import completed fem.py (to complete)
import fem
import plot


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


class TUTangent(TU):

    FUNCTIONS = ["new_file",
                 "new_material",
                 "open_mesh",
                 "create_model",
                 "calculate_vectorfunction",
                 "compute_tangent",
                 "compute_Knum",
                 "compare_results"]

    FILENAMES = ["triangle-quad11.msh",
                 "triangle-quad44.msh"]
    # FILENAMES = ["triangle-quad11.msh"]

    MATERIALS = [elasticity.StVenantKirchhoffElasticity,
                 elasticity.NeoHookeanElasticity]

    # MATERIALS = [elasticity.StVenantKirchhoffElasticity]

    def __init__(self):
        super(TUTangent, self).__init__()

    def __setUp(self):
        number_files = len(TUTangent.FILENAMES)
        for i in range(number_files):
            TUTangent.FILENAMES[i] = TU.FOLDER_MSH + TUTangent.FILENAMES[i]
        number_materials = len(TUTangent.MATERIALS)

        TUTangent.FUNCTIONS *= number_materials * number_files
        self.set_DEN()
        self.n_mat = 0
        self.n_fil = 0

    def new_file(self):
        """
        Here we change the filename to test
        """
        self.filename = TUTangent.FILENAMES[self.n_fil]
        self.n_fil += 1

    def new_material(self):
        if self.n_mat == 0:
            E = 210.e9
            nu = 0.3
        else:
            E = 10.e6
            nu = 0.45
        self.material = TUTangent.MATERIALS[self.n_mat](E, nu)
        if self.n_fil >= len(TUTangent.FILENAMES):
            self.n_fil = 0
            self.n_mat += 1

    def open_mesh(self):
        meshfile = open(self.filename, 'r')
        self.mesh = gmsh.gmshInput_mesh(meshfile)
        meshfile.close()

    def create_model(self):
        self.FE_model = fem.FEModel(self.mesh, self.material)

    def calculate_vectorfunction(self):
        """
        Calculate the vector deplacement vector U using the function "fct"
        """
        def fct(x):
            fx = x[0] * x[1]
            fy = x[0] * x[1]
            return [fx, fy]

        nNodes = self.mesh.getNNodes()
        dim = self.mesh.getDim()
        self.U = np.zeros((nNodes, dim))
        for i in range(nNodes):
            Nodei = self.mesh.getNode(i)
            Xi = Nodei.getX()
            self.U[i] = fct(Xi)

    def compute_tangent(self):
        """
        get gradient and stress tensor fields
        """
        self.K = self.FE_model.computeTangent(self.U)

    def compute_Knum(self):
        DU = 1.e-5  # perturbation on the displacements U
        # Rp = np.zeros((mesh.nNodes(), 2))
        # Rm = np.zeros((mesh.nNodes(), 2))

        nNodes = self.mesh.getNNodes()
        dim = self.mesh.getDim()
        self.Knum = np.zeros((nNodes, dim, nNodes, dim))
        Ud = np.copy(self.U)

        for n in range(nNodes):
            for j in range(dim):
                Ud[n, j] += DU / 2.
                # R(U_nj + DU/2), residual "plus"
                Rp = self.FE_model.computeResidual(Ud)

                Ud[n, j] -= DU
                # R(U_nj - DU/2), residual "minus"
                Rm = self.FE_model.computeResidual(Ud)

                # come back to before perturbation for next iterations
                Ud[n, j] += DU / 2

                self.Knum[n, j] = (Rp - Rm) / DU
                # for m in range(nNodes):
                #     for i in range(dim):
                #         # self.Knum[n, j, m, i] = Ktemp[m, i]
                #         self.Knum[m, i, n, j] = Ktemp[m, i]

    def compare_results(self):
        Ktest = np.absolute(np.copy(self.K))
        Ktest[Ktest < 1] = 1
        diffK = np.fabs(self.K - self.Knum) / Ktest
        KerrMax = np.max(diffK)

        TOLMAX = 1.e-2  # absolute tolerance on error
        if KerrMax > TOLMAX:
            msg = "Code NOT VALIDATED: [KerrMax/TOLMAX] = [%.2e/%.2e]" % (
                KerrMax, TOLMAX)
            raise Exception(msg)


if __name__ == "__main__":
    test = TUTangent()
    result = test.run()
