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


#
#--- read mesh from Gmsh file
#


msh_folder = "../msh/"
basename = "triangle-quad11.msh"
# basename = "triangle-quad44.msh"
inputfnam = msh_folder + basename
meshfile = open(inputfnam, 'r')
mesh = gmsh.gmshInput_mesh(meshfile)
meshfile.close()

print("Read mesh '%s' with %d nodes and %d elements" %
      (inputfnam[:-4], mesh.nNodes(), mesh.nElements()))

#
#--- constitutive model and parameters (test both !)
#
#E = 210.e9
#nu = 0.3
#mat_model = elasticity.StVenantKirchhoffElasticity(E,nu)
E = 10.e6
nu = 0.45
mat_model = elasticity.NeoHookeanElasticity(E, nu)
print("Chose material model '%s' with E = %.2e Pa and nu = %.2f" %
      (mat_model.__class__.__name__, E, nu))

#
#--- create FE model
#

FE_model = fem.FEModel(mesh, mat_model)
nElements = FE_model.getNElements()
nNodes = mesh.getNNodes()
dim = mesh.dim

#
#--- define vector function
#


def fct(x):
    fx = x[0] * x[1]
    fy = x[0] * x[1]
    return [fx, fy]


#
#--- create nodal array
#
U = np.zeros((nNodes, dim))  # initialise array of zeros,
# of shape nNodes x 2: [U1_x,U1_y,...,UnNodes_x,UnNodes_y]
for i in range(mesh.nNodes()):
    Nodei = mesh.getNode(i)
    Xi = Nodei.getX()
    U[i] = fct(Xi)

#
#--- compute tangent matrix
#
# Nnodes = mesh.nNodes()
# dim = 2
# K = np.zeros((Nnodes, dim, Nnodes, dim))
K = FE_model.computeTangent(U)
# print("K = ")
# print(K)

# =============================================================================================
# Compute tangent by perturbation (finite central difference) and compare to analytical
# =============================================================================================
# K_minj  = dR_mi / dU_nj
# hence Knum_minj = ( R(U_nj + DU/2)_mi - R(U_nj - DU/2)_mi )/DU
# no need to modify this part, it is for help

print("\nComputing error on tangent matrix...\n")

# error printed in logfile or in console
# whether or not to print in console element wise error between K and Knum (line 128)
printerror = False
logfilename = inputfnam[:-4] + "-test_tangent.txt"
logfile = open(logfilename, 'w')
logfile.write('---INFO---: m i n j  K[m, i, n, j]   Knum[m, i]      Kerr\n')

Ud = np.array(U)
DU = 1.e-5  # perturbation on the displacements U
TOLMAX = 1.e-2  # absolute tolerance on error
PRECISION = 1e-16  # zero precision
# Rp = np.zeros((mesh.nNodes(), 2))
# Rm = np.zeros((mesh.nNodes(), 2))
KerrMax = 1e-10


for n in range(nNodes):
    for j in range(dim):
        Ud[n, j] += DU / 2.
        Rp = FE_model.computeResidual(Ud)  # R(U_nj + DU/2), residual "plus"

        Ud[n, j] -= DU
        Rm = FE_model.computeResidual(Ud)  # R(U_nj - DU/2), residual "minus"

        Ud[n, j] += DU / 2  # come back to before perturbation for next iterations

        Knum = (Rp - Rm) / DU
        # print("Knum = " + str(Knum))
        for m in range(nNodes):
            for i in range(dim):
                # Knum = (Rp[m, i] - Rm[m, i]) / DU
                # print("Knum = " + str(Knum.shape))
                Ktest = max(1, np.fabs(K[m, i, n, j]))
                # print("np.fabs = " + str(np.fabs(K[m, i, n, j] - Knum[m, i])))
                Kerr = np.fabs(K[m, i, n, j] - Knum[m, i]) / Ktest
                KerrMax = max(KerrMax, Kerr)
                # if (Kerr > KerrMax):
                #     KerrMax = Kerr
                if (printerror) and (Kerr > TOLMAX) and (np.fabs(Knum[m, i]) > PRECISION):
                    print('%d %d %d %d %13.6e %13.6e %13.6e' %
                          (m, i, n, j, K[m][i][n][j], Knum[m, i], Kerr))
                if (Kerr > TOLMAX):  # and (np.fabs(Knum)>PRECISION) :
                    logfile.write('-CRITICAL-: %d %d %d %d %13.6e %13.6e %13.6e\n' %
                                  (m, i, n, j, K[m, i, n, j], Knum[m, i], Kerr))
                else:
                    logfile.write('---INFO---: %d %d %d %d %13.6e %13.6e %13.6e\n' %
                                  (m, i, n, j, K[m, i, n, j], Knum[m, i], Kerr))
        # print("flag = " + str(flag))
        # print("Knum = ")
        # print(Knum)
        # stop()


# Print result
print('Maximal error = %13.6e' % KerrMax)
print("Tolerance = %13.6e" % TOLMAX)
line = "***************************************************************\n\
**************** Maximal error =  %13.6e ***************\n" % KerrMax
logfile.write(line)

if KerrMax > TOLMAX:
    print("!!!!!!!!!!!! Code NOT VALIDATED: KerrMax = %.2e !!!!!!!!!!!!" % KerrMax)
    logfile.write("!!!!!!!!!!!! Code NOT VALIDATED: KerrMax = %.2e > TOLMAX = %.0e!!!!!!!!!!!!" % (
        KerrMax, TOLMAX))
else:
    print("************ Code VALIDATED *************")
    line = "*********************** Code VALIDATED ************************\n"
    line += "***************************************************************"
    logfile.write(line)
logfile.close()
