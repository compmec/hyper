# -*- coding: utf-8 -*-
#
# Script to test consitutive models implemented in
#


"""
MNMNL Homework
Module 4: stress tensor field on a mesh

This script tests your implementation of the elasticity.py module with numerical
differentiation of the elastic potential for the Piola-Kirchhoff II stress
tensor S and of the stress for the material tangent operator M.

You will test two types of constitutive models:
    - StVenantKirchhoff compressible material
    - NeoHookean compressible material

You will test two types of meshes (with two refinements):
    - triangle-tri(14,56).msh with 14 and 56 triangular elements respectively
    - meshfile-quad(11,44).msh with 11 and 44 quadrangular elements respectively

You can also compare your stress tensor field with the anaylitical solution you
can derive from the exercise.

WARNING: THERE IS NOTHING TO COMPLETE HERE ! YOU MUST COMPLETE M4__elasticity_template.py
Simply comment/uncomment meshfiles (lines 51-54) and constitutive models (lines 66-67)
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
    import tensor
    import elasticity
except ModuleNotFoundError as e:
    error = str(e)
    notimportedfilename = error.split("'")[1]
    print("Import file not found: " + str(notimportedfilename) + ".py")
    sys.exit()

try:
    TU.include_path("msh")
except ModuleNotFoundError:
    print("Could not include the folder with the mesh examples")


class TUElasticity(TU):

    FUNCTIONS = ["create_nodalarray", "compute_stresstensorfield",
                 "write_fields", "compute_stress", "compute_materialtangent"]

    def __init__(self, inputfilename, material):
        self.filename = inputfilename
        self.material = material
        super(TUElasticity, self).__init__()

    def __setUp(self):
        meshfile = open(self.filename, 'r')
        self.mesh = gmsh.gmshInput_mesh(meshfile)
        meshfile.close()

        print("Read mesh with %d nodes and %d elements" %
              (self.mesh.nNodes(), self.mesh.nElements()))

        # Create FE Model
        self.FE_model = fem.FEModel(self.mesh, self.material)

    def create_nodalarray(self):
        def fct(x):
            f0 = x[0] * x[1]
            f1 = x[0] * x[1]
            return [f0, f1]

        mesh = self.mesh
        self.U = np.zeros((mesh.nNodes(), 2))
        for i in range(mesh.nNodes()):
            self.U[i] = fct([mesh.getNode(i).getX(0), mesh.getNode(i).getX(1)])

    def compute_stresstensorfield(self):
        #
        #--- compute stress tensor field
        #

        self.FE_model.computeStress(self.U)
        self.F = []
        self.E = []
        self.P = []
        nElements = self.FE_model.nElements()
        for n in range(nElements):  # tensor at each element of the mesh
            new_F = self.FE_model.getDeformationGradient(n)
            new_E = self.FE_model.getStrainGreenLagrange(n)
            new_P = self.FE_model.getStressPK1(n)
            self.F.append(new_F.flatten())
            self.E.append(new_E.flatten())
            self.P.append(new_P.flatten())

    def write_fields(self):
        #
        #--- write fields
        #
        # filename of the new meshfile
        outputfnam = self.filename[:-4] + '-val.msh'
        outfile = open(outputfnam, 'w')
        gmsh.gmshOutput_mesh(outfile, self.mesh)
        # creates a view in gmsh, vector field
        gmsh.gmshOutput_nodal(outfile, "displacement U", self.U, 0, 0.0)
        # creates a view in gmsh, tensor field
        gmsh.gmshOutput_element(
            outfile, "deformation gradient F", self.F, 0, 0.0)
        # creates a view in gmsh, tensor field
        gmsh.gmshOutput_element(
            outfile, "engineering stress P", self.P, 0, 0.0)
        outfile.close()

    def compute_stress(self):

        # =============================================================================================
        #%% Compute error on stress by finite (central) difference of energy potential
        # =============================================================================================
        d = self.FE_model.getDim()
        DC = 1.e-5
        TOLMAX = 1.e-3
        PRECISION = 1.e-16
        printerror = True
        e = 1  # element
        SerrMax = 1e-30
        nElements = self.FE_model.nElements()
        for e in range(nElements):
            Ed = self.E[e].reshape(d, d)  # Lagrange deformation
            Cd = 2 * Ed + tensor.I(d)  # right Cauchy-Green tensor
            S = self.FE_model.material.stress(Cd)
            for i in range(d):
                for j in range(d):

                    Cd[i, j] += DC / 2.
                    # phi(C_ij + DC/2), potential plus
                    phip = self.FE_model.material.potential(Cd)

                    Cd[i, j] -= DC
                    phim = self.FE_model.material.potential(Cd)

                    Cd[i, j] += DC / 2.  # come back to before perturbation

                    # S_ij = 2*(phi(C_ij+DC/2) - phi(C_ij-DC/2))/DC
                    Snum = 2 * (phip - phim) / DC

                    if (np.fabs(S[i, j]) > 1.0):
                        Stest = np.fabs(S[i, j])
                    else:
                        Stest = 1.0
                    Serr = np.fabs(S[i, j] - Snum) / Stest
                    if (Serr > SerrMax):
                        SerrMax = Serr
                    if (printerror) and (Serr > TOLMAX) and (np.fabs(Snum) > PRECISION):
                        print('%d %d %13.6e %13.6e %13.6e' %
                              (i, j, S[i, j], Snum, Serr))

        # print('Maximal error = %13.6e' % SerrMax)
        if SerrMax > TOLMAX:
            msg = "Stress computation not valided: SerrMax = %.2e" % SerrMax
            raise Exception(msg)

    def compute_materialtangent(self):
        # =============================================================================================
        #%% Compute error on material tangent operator by finite (central) difference of stress
        # =============================================================================================
        d = self.FE_model.dim
        DC = 1.e-5
        TOLMAX = 1.e-3
        PRECISION = 1.e-16
        printerror = False
        e = 1  # element
        KerrMax = 1e-30
        for e in range(self.FE_model.nElements()):
            Ed = self.E[e].reshape(d, d)  # Lagrange deformation
            Cd = 2 * Ed + tensor.I(d)  # right Cauchy-Green tensor
            Sp = tensor.tensor(d)  # S(C+DC/2) Piola Kirchhoff II tensor "plus"
            # S(C-DC/2) Piola Kirchhoff II tensor "minus"
            Sm = tensor.tensor(d)
            S = self.FE_model.material.stress(Cd)
            # Lagrangian tangent operator 2*dS/dC
            M = self.FE_model.material.stiffness(Cd)

            for k in range(d):
                for l in range(d):
                    Cd[k, l] += DC / 2.
                    Sp = self.FE_model.material.stress(Cd)

                    Cd[k, l] -= DC
                    Sm = self.FE_model.material.stress(Cd)

                    Cd[k, l] += DC / 2.  # come back to before perturbation
                    for i in range(d):
                        for j in range(d):
                            Knum = (Sp[i, j] - Sm[i, j] + Sp[j, i] -
                                    Sm[j, i]) / DC  # 2*(Sp-Sm + Sp-Sm)/2/DC
                            # K_ijkl = ( S_ij(C_kl+DC/2) - S_ji(C_kl-DC/2) ) SYM
                            if (np.fabs(M[i, j, k, l]) > 1.0):
                                Ktest = np.fabs(M[i, j, k, l])
                            else:
                                Ktest = 1.0
                            Kerr = np.fabs(M[i, j, k, l] - Knum) / Ktest
                            if (Kerr > KerrMax):
                                KerrMax = Kerr
                            if (printerror) and (Kerr > TOLMAX) and (np.fabs(Knum) > PRECISION):
                                print('%d %d %d %d %13.6e %13.6e %13.6e' %
                                      (i, j, k, l, M[i, j, k, l], Knum, Kerr))
        # print('Maximal error = %13.6e' % KerrMax)
        if KerrMax > TOLMAX:
            msg = "Material tangent not valided: KerrMax = %.2e" % KerrMax
            raise Exception(msg)


if __name__ == "__main__":

    msh_folder = "../msh/"
    files = ["triangle-tri14.msh", "triangle-tri56.msh",
             "triangle-quad11.msh", "triangle-quad44.msh"]
    materials = [elasticity.StVenantKirchhoffElasticity,
                 elasticity.NeoHookeanElasticity]

    for i in range(len(materials)):
        for j in range(len(files)):
            print("##############################")
            print("#   Test Material Number " + str(i + 1) + "   #")
            print("#   Test File Number " + str(j + 1) + "       #")
            print("##############################")

            E = 210.e9  # Pa
            nu = 0.3
            mat = materials[i](E, nu)
            filename = msh_folder + files[j]
            test = TUElasticity(filename, mat)
            result = test.run()
            if result == test.FAILURE:
                break
