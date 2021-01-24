# -*- coding: utf-8 -*-
#
# Tensile test on a membrane with hole, assuming plane strain conditions
#

"""
MNMNL Homework
Module 7: non linear solver and 2D plane strain problem

Complete code below where indicated by
# ============
# TODO:
# ...
# ============

In this script, you must complete the non linear Newton-Raphson solver to compute
the results for a pressure test on a compressible hyperelastic thick cylinder,
in plane strain conditions (2D).
There is no external forces: the cylinder is loaded by means of prescribed radial
displacements on the interior surface. (hence R = Tint and K = Kint as in previous homeworks)

You must test BOTH material models:
    - St Venant Kirchhoff
    - NeoHookean

You must COMPULSORILY test both triangular meshes:
    - cynlinder-tri_coarse.msh
    - cynlinder-tri_refine.msh

You can also test quandragular meshes:
    - cynlinder-quad_coarse.msh
    - cynlinder-quad_refine.msh

In the post-processing section, there is nothing to complete but you can modify
it as you wish.
It provides 2 figures, save in .png files, showing the loading curve of the test
and the evolution of the displacement along the bottom line, to give a comparison
with the linear elastic result:
    u = u_r * e_r = (a*r + b/r) * e_r
with (a,b) constants which depend on the loading and the material properties.

"""
import os
import sys
sys.path.append("../src/")
sys.path.append("../geo/")
sys.path.append("../msh/")
# os.chdir('D:\\Documents\\Hyper\\Homeworks\\HW7')
import numpy as np
import matplotlib.pyplot as plt
import gmsh2 as gmsh
import fem
import elasticity
from scipy.linalg import logm
from datetime import datetime
import errno
import plot


def mkdir_p(path):
    """ make directory at 'path' (full or relative) if not existing """
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def stop():
    raise Exception("stop")


def cut(M, bcDofs, freeDofs):
    """
    Function to cut a numpy array, like:
        If M is a vector like, we cut it in M1 and M2 like:
            M1 = M[bcDofs]
            M2 = M[freeDofs]
        If M is a matrix, we cut like:
            M = [ M11  M12 ]
                [ M21  M22 ]
            M11 = M[bcDofs, bcDofs]
            M12 = M[bcDofs, freeDofs]
            M21 = M[freeDofs, bcDofs]
            M22 = M[freeDofs, freeDofs]

    """
    if len(M.shape) == 1:
        # F1 = np.delete(M, freeDofs)
        # F2 = np.delete(M, bcDofs)
        F1 = M[bcDofs]
        F2 = M[freeDofs]
        return (F1, F2)
    else:
        M1 = np.delete(M, freeDofs, axis=0)
        M2 = np.delete(M, bcDofs, axis=0)
        M11 = np.delete(M1, freeDofs, axis=1)
        M21 = np.delete(M2, freeDofs, axis=1)
        M12 = np.transpose(M21)
        M22 = np.delete(M2, bcDofs, axis=1)
        # M11 = M[bcDofs, bcDofs]
        # M12 = M[bcDofs, freeDofs]
        # M21 = M[freeDofs, freeDofs]
        # M22 = M[freeDofs, freeDofs]
        return (M11, M12, M21, M22)


def verify(U, Uref, bcDofs, TOLERANCE=1e-9):
    """
    This function verify if the BC conditions always stays the same.
    As we know that U[bcDofs] doesn't change, during the operations of U, we verify if it has changed using the Uref.
    Uref has the good conditions at bcDofs positions, and 0 in the rest
    """
    if norme((U[bcDofs] - Uref[bcDofs])) > TOLERANCE:
        print("norm(U[freeDofs]) = " + str(norme(U[freeDofs])))
        print("U[bcDofs] = " + str(U[bcDofs]))
        print("Uref[bcDofs] = " + str(Uref[bcDofs]))
        raise Exception("U[bcDofs] != Uref[bcDofs]")


print('============================================')
print('=                                          =')
print('=      2D plane strain thick cylinder      =')
print('=                                          =')
print('============================================')


VERBOSE = True

###############################################################################
#                PROBLEM DEFINITION                                           #
###############################################################################

# DEFAULT UNITS: m, Pa, N

#################################
#            MATERIAL           #
#################################
material = "steel"
matmodel = 'StVenant'
# material = "rubber"
# matmodel = 'NeoHookean'
if material == "steel":
    # steel
    Eyoung = 210e5  # Young's modulus (Pa)
    nu = 0.3  # Poisson's ratio
else:
    # filled rubber (very stiff):
    Eyoung = 10e6  # Young's modulus (Pa)
    nu = 0.45  # Poisson's ratio

if VERBOSE:
    print("############################")
    print("#         MATERIAL         #")
    print("############################")
    print("    Material: " + str(material))
    print("           E = %.2e" % Eyoung)
    print("           v = %.2e" % nu)
    print("       model = " + str(matmodel))
    print("")


#################################
#      ITERATION PARAMETERS     #
#################################
nSteps = 10  # number of steps in the incremental loading loop
tStep = 1. / nSteps  # size of the steps
# rdisp = 0.00001  # constant radial displacement on interior boundary (StVenant)
rdisp = 0.0001  # constant radial displacement on interior boundary (StVenant)
# rdisp = 0.01 #constant radial displacement on interior boundary (NeoHookean)
if VERBOSE:
    print("############################")
    print("#   ITERATION PARAMETERS   #")
    print("############################")
    print("     Number of steps = " + str(nSteps))
    print("             delta T = " + str(tStep))
    print(" radial displacement = %.2e" % rdisp)
    print("")


#################################
#        SOLVER PARAMETERS      #
#################################
PRECISION = 1e-15  # machine "zero"
TOLERANCE = 1e-6
ITERMAX = 10
norme = np.linalg.norm
if VERBOSE:
    print("############################")
    print("#     SOLVER PARAMETERS    #")
    print("############################")
    print("    PRECISION = %.2f" % PRECISION)
    print("    TOLERANCE = %.2f" % TOLERANCE)
    print("      ITERMAX = %.2f" % ITERMAX)
    print("")


#################################
#       FILES AND FOLDERS       #
#################################

geo_folder = "../geo/"
msh_folder = "../msh/"

basemesh = "cylinder-custom.msh"
# basemesh = "cylinder-tri_coarse.msh"
# basemesh = "cylinder-tri_refined.msh"
# basemesh = "cylinder-quad_coarse.msh"
# basemesh = "cylinder-quad_refined.msh"

basegeo = 'cylinder.geo'
if VERBOSE:
    print("############################")
    print("#     FILES AND FOLDERS    #")
    print("############################")
    print("    mesh = " + str(basemesh))
    print("    geom = " + str(basegeo))


#################################
#    READ FILES INFORMATIONS    #
#################################

material_folder = msh_folder + matmodel
dirout = os.path.join(msh_folder, matmodel)

with open(geo_folder + basegeo) as geomfile:
    lines = geomfile.readlines()
    Ri = float(lines[0].strip().strip(";").split(
        " = ")[1])  # interior radius (m)
    Re = float(lines[1].strip().strip(";").split(
        " = ")[1])  # exterior radius (m)

if VERBOSE:
    print("    Geometry cylinder:")
    print("        Internal radius Ri = %.2e m" % Ri)
    print("        External radius Re = %.2e m" % Re)


inputfile = msh_folder + basemesh
meshfile = open(inputfile, 'r')
mesh = gmsh.gmshInput_mesh(meshfile, d=2)
meshfile.close()
dim = mesh.dim  # dimension of the problem
nNodes = mesh.nNodes()  # number of nodes in the mesh
nElements = mesh.nElements()  # number of elements in the mesh
if VERBOSE:
    print("    Mesh read:")
    print("        File: %s" % inputfile)
    print("        Nodes: %d" % nNodes)
    print("        Elements: %d" % nElements)
    print("")


#################################
#          GSMH OUTPUT          #
#################################

mesh_name = basemesh.split('-')[0]
# 'tri_coarse' = element type + mesh refinement
mesh_type = basemesh.split('-')[1][:-4]
# subdirectory for output files
new_folder_name = mesh_type

dirout = os.path.join(msh_folder + matmodel, new_folder_name)
mkdir_p(dirout)  # create directory for output files if not already existing
outputfile = mesh_name + '-results_Nsteps-' + str(int(nSteps)) + ".msh"
outputfile = os.path.join(dirout, outputfile)
gmsh_out = open(outputfile, 'w')
gmsh.gmshOutput_mesh(gmsh_out, mesh)

plot.plot_mesh(mesh, label="Before", color="b")


#################################
#         PREPROSSESING         #
#################################

#--- boundary conditions
#
# On the left side of the quarter-cylinder (X=0), the displacements in the x-direction are fixed (u_x =0)
# At the bottom of the quarter-cylinder (Y=0), the displacements in the y-direction are fixed (u_y=0)
# We apply a radial displacement "disp" on the interior of the cylinder, in the r-direction (u_r=rdisp)

nNodes = mesh.nNodes()
nDofs = dim * nNodes

# =============================================================================
# TODO: DONE and verified with mesh file
# list degrees of freedom (d.o.f) in x-direction linked to left side and
# d.o.f in y-direction linked to bottom side in the list 'fixDofs'
# Hint: for loop on all the nodes and check coordinates with mesh.getNode(i).getX()
Xmin = np.min([mesh.getNode(inode).getX(0) for inode in range(nNodes)])
# Maillage pas forcément centré en (0,0)
Ymin = np.min([mesh.getNode(inode).getX(1) for inode in range(nNodes)])


# =============================================================================
leftNodes = []
bottomNodes = []
dispNodes = []


normals = []

for i in range(nNodes):
    Nodei = mesh.getNode(i)
    X, Y = Nodei.getX()[:2]
    if np.abs(X) < PRECISION and np.abs(Y) > Ri - PRECISION:
        leftNodes.append(i)
    if np.abs(Y) < PRECISION and np.abs(X) > Ri - PRECISION:
        bottomNodes.append(i)

for i in range(nNodes):
    Nodei = mesh.getNode(i)
    X, Y = Nodei.getX()[:2]
    if np.abs(X**2 + Y**2 - Ri**2) < PRECISION:  # X^2 + Y^2 = Ri^2
        dispNodes.append(i)
        normals.append([X, Y])


left_fixDofs = []  #
bottom_fixDofs = []
dispDofs = []  # initialise empty list
for i in leftNodes:
    left_fixDofs.append(i * dim)
for i in bottomNodes:
    bottom_fixDofs.append(i * dim + 1)
for i in dispNodes:
    dispDofs.append(i * dim)
    dispDofs.append(i * dim + 1)
fixDofs = left_fixDofs + bottom_fixDofs
fixDofs = list(np.sort(fixDofs))
dispDofs = list(np.sort(dispDofs))


normals = np.array(normals) / Ri
# applied to every node according normal to the surface in that node
dispInt = rdisp * normals
# from 2-entry table to flattened array [u1,v1,u2,v2,...,uN,vN]
disp = dispInt.flatten()
# d.o.f. subjected to a boundary condition:
bcDofs = fixDofs + dispDofs  # the addition of 2 lists concatenates them
# d.o.f. free of any constraint (where the displacement field is to be determined):
bcDofs = list(dict.fromkeys(bcDofs))
bcDofs = list(np.sort(bcDofs))

listDofs = np.arange(nDofs)
# remove all entries in bcDofs from listDofs
freeDofs = np.delete(listDofs, bcDofs)
freeDofs = list(np.sort(freeDofs))
# display fixed degrees of freedom
if VERBOSE:
    print('############################')
    print('#      PRE PROCESSING      #')
    print('############################')
    print("    Mesh information:")
    print("        Total nodes: %d" % nNodes)
    print("        Total elements: %d" % nElements)
    print("        Dimention: %d" % dim)
    print("    Number of nodes:")
    print("        Left: %d" % len(leftNodes))
    print("        Center: %d" % len(dispNodes))
    print("        Bottom: %d" % len(bottomNodes))
    print("    Degrees of freedom:")
    print("        Total: %d" % nDofs)
    print("        BC Left: %d" % len(bottom_fixDofs))
    print("        BC Center: %d" % len(dispDofs))
    print("        BC Bottom: %d" % len(bottom_fixDofs))
    print("        BC Total: %d" % len(bcDofs))
    print("        Free: %d" % len(freeDofs))
    print("")

    # print('    Fixed degrees of freedom: [' +
    #   ', '.join(['{:d}'.format(fixdof) for fixdof in fixDofs]) + ']')
# bcDofs = list(set(bcDofs)) #If there are duplicates in bcDofs

#
# stop()

#
#--- create material model instance
#
if matmodel == "StVenant":
    material = elasticity.StVenantKirchhoffElasticity(Eyoung, nu)
elif matmodel == "NeoHookean":
    material = elasticity.NeoHookeanElasticity(Eyoung, nu)
else:
    raise(KeyError, "matmodel must be one of 'StVenant', 'NeoHookean'")

#
#--- create FE model instance
#
model = fem.FEModel(mesh, material)

#################################
#       INITIALISE ARRAYS       #
#################################
# displacement array
# flattened vector of the 2-entry array U = [u1 v1 u2 v2 ... uN vN]
U_notflat = np.zeros((nNodes, dim))
U = np.zeros(nDofs)
U1 = np.zeros(len(bcDofs))
U2 = np.zeros(len(freeDofs))
Uref = np.zeros(nDofs)
# residual array
Fint_notflat = np.zeros((nNodes, dim))  # unknown: initialisation (R=Tint),
Fext_notflat = np.zeros((nNodes, dim))  # unknown: initialisation (R=Tint),
Fint = np.zeros(nDofs)  # unknown: initialisation (R=Tint),
F1 = np.zeros(len(bcDofs))
F2 = np.zeros(len(freeDofs))
Fext = np.zeros(nDofs)  # unknown: initialisation (R=Tint),
# flattened vector of the 2-entry array U = [Rx1 Ry1 Rx2 Ry2 ... RxN RyN]
# tangent matrix
# flattened matrix of the 4-entry array K = [K11.., K12.., K21.., K22.., ..., KN1.., KN2..]
K = np.zeros((nDofs, nDofs))


# R_not_flat = model.computeResidual(U.reshape(-1,2),R.reshape(-1,2))
# R = R_not_flat.flatten()
# print('WESH Norme de R : ', np.linalg.norm(R))


#################################
#           SIMULATION          #
#################################
if VERBOSE:
    print('############################')
    print('#        SIMULATION        #')
    print('############################')
start_time = datetime.now()  # measure time for the simulation
time = 0.0  # final time is nSteps x tStep = nSteps x (1/nSteps) = 1
itermoy = 0.  # average number of iterations per time step
# prescribed displacement per step on the right hand side of the membrane
displacements = []
reactions = []  # corresponding reaction per step

for iStep in range(nSteps):  # loop on time steps --> incremental loading
    time += tStep
    if VERBOSE:
        print("--> Step [%d/%d]" % (iStep + 1, nSteps))
        print("    time = %.2f/1.00" % time)
        print("    disp = %.2e" % (rdisp * time))

    # =========================================================================
    #
    # We know some displacements, but not all. Let's say that we only know
    # the vector displacement U1 and the forces F2 that are applied in the
    # other's nodes. So, if we have the system:
    #
    #     [ M ] * { U } = { F }
    #
    # We can solve as
    #
    #     [ [M11]   [M12] ]   { {U1} }   { {F1} }
    #     [               ] * {      } = {      }
    #     [ [M21]   [M22] ]   { {U2} }   { {F2} }
    #
    # So, we can write the unknown solution of U2 like:
    #
    #   U2 = M22_inv * ( F2 - M21 * U1 )
    #
    # and after that, the reaction forces:
    #
    #   F1 = M11 * U1 + M12 * U2
    #
    # So, we know:
    #     * {U1} - the central displacement
    #     * {F2} - the forces aux noeuds that are not in Boundary Conditions
    # And the steps to calculate are:
    #     * Calculate all matrix [M] using {U} = {0} in the model
    #     * Divide [M] into 4 matrix: [M11], [M12], [M21], [M22]
    #     * Calculate {U2} using {U1} and {F2}
    #     * Mix {U1} and {U2} to create the displacement vector {U}
    #     * Calculate the internal force {Fint} using {U} in the model
    #     * Calculate {F1} = {Fext} in the BC that makes the displacement {U1}
    #     * Mix {F1} and {F2} to create the external force {Fext}
    #     * Compute the Residu {R} = {Fext} - {Fint}
    #     * We cut {R} into {R1} and {R2}
    #     * We compare the value of the norm {R2} (BC region):
    #         - If it's less than a tolerance, that's good
    #         - If it's bigger than the tolerance, we do:
    #         - Calculate [M] again, with the new value of {U}
    #         - As the values of {U1} doest change, we compute
    #
    #             dU2 = M22_inv * R2
    #
    #         - Sum {dU2} in {U2} and put it in the {U}
    #         - Compute the internal force again {Fint} using the model
    #         - Compute {dF1} using the value of {dU2} and put it in {Fext}
    #         - Compute {R} = {Fext} - {Fint} and cut it in {R1} and {R2}
    #         - Compare the value of the norm {R2}, return if the tolerance
    #             is not good enough.
    #
    # =========================================================================

    K_notflat = model.computeInternalStiffness(U_notflat)
    K = K_notflat.reshape((nDofs, nDofs))
    (M11, M12, M21, M22) = cut(K, bcDofs, freeDofs)

    Fint_notflat = model.computeInternalForces(U_notflat)
    Fint = Fint_notflat.flatten()
    (F1, F2) = cut(Fint, bcDofs, freeDofs)

    for j, idof in enumerate(bcDofs):
        if idof in dispDofs:
            i = dispDofs.index(idof)
            U1[j] = disp[i] * time
            Uref[idof] = disp[i] * time

    U2 = np.linalg.solve(M22, F2 - M21 @ U1)
    U[bcDofs] = U1[:]
    U[freeDofs] = U2[:]
    U_notflat = U.reshape((nNodes, dim))

    verify(U, Uref, bcDofs)  # Verify the good BC conditions
    Fext[bcDofs] = M11 @ U1 + M12 @ U2

    R = Fext - Fint  # Residu
    (R1, R2) = cut(R, bcDofs, freeDofs)

    # =========================================================================

    norm = np.linalg.norm(R)
    test = TOLERANCE * norm
    if (test < PRECISION):
        test = PRECISION
    iteration = 0
    if VERBOSE:
        print('    *** Iteration %02d: residual norm =%.3e ***' %
              (iteration, norm))

    while (norm > test):  # Newton-Raphson iterations

        # compute correction
        # (U1, U2) = cut(U, bcDofs, freeDofs)
        # U_notflat = U.reshape((nNodes, dim))

        K_notflat = model.computeInternalStiffness(U_notflat)
        K = K_notflat.reshape((nDofs, nDofs))
        (M11, M12, M21, M22) = cut(K, bcDofs, freeDofs)

        dU2 = np.linalg.solve(M22, R2)
        U[freeDofs] += dU2[:]
        U_notflat = U.reshape((nNodes, dim))

        Fint_notflat = model.computeInternalForces(U_notflat)
        Fint = Fint_notflat.flatten()

        Fext[bcDofs] = M11 @ U1 + M12 @ U2

        R = Fext - Fint
        R1, R2 = cut(R, bcDofs, freeDofs)

        verify(U, Uref, bcDofs)  # Verify the good BC conditions

        normdU2 = np.linalg.norm(dU2)
        if (normdU2 < PRECISION) and VERBOSE and iteration == ITERMAX:
            print("    ERROR: no convergence in Newton loop, dQ stagnation")

        norm = np.linalg.norm(R2)
        iteration = iteration + 1
        if VERBOSE:
            print('    *** Iteration %02d: residual norm =%.3e ***' %
                  (iteration, norm))

        # =====================================================================
        if iteration == ITERMAX:
            print("    *** ITERMAX reached for step: %02d" % iStep)
            print("    ***            Residual norm: %.3e " % norm)
            break
        # =====================================================================
    # end of Newton-Raphson iterations
    itermoy += iteration
    if norm <= test and VERBOSE:
        print("    COMPLETE: norm < test")
        print("        norm = %.5e" % norm)
        print("        test = %.5e" % test)

    nodalForce = []
    uint = []
    for i in dispNodes:
        u = np.array(U.reshape(nNodes, dim)[i, :]).flatten()
        X = np.array([mesh.getNode(i).getX(0), mesh.getNode(i).getX(1)])
        x = X + u
        e_r = x / np.linalg.norm(x)
        nodalForce.append(R.reshape(nNodes, dim)[i, :].dot(e_r))
        uint.append(u.dot(e_r))
        # assert np.allclose(rdisp*time,uint[-1])
    reactions.append(nodalForce)  # (N)
    displacements.append(uint)  # (m)
    #
    #--- get strain and stress tensor fields for current time step
    #
    E = []
    P = []
    hencky = []
    euler = []
    sigma = []
    for n in range(model.nElements()):
        E.append(model.getStrainGreenLagrange(n).flatten())
        P.append(model.getStressPK1(n).flatten())
        hencky.append(model.getStrainHencky(n).flatten())
        euler.append(model.getStrainEulerAlmansi(n).flatten())
        sigma.append(model.getStressCauchy(n).flatten())
    #
    #--- gmsh output for current time step
    #
    gmsh.gmshOutput_nodal(gmsh_out, "Residual",
                          R.reshape(nNodes, dim), iStep, time)
    gmsh.gmshOutput_nodal(gmsh_out, "Displacements",
                          U.reshape((nNodes, dim)), iStep, time)
    gmsh.gmshOutput_element(gmsh_out, "Green-Lagrange strain", E, iStep, time)
    gmsh.gmshOutput_element(
        gmsh_out, "Piola Kirchhoff I stress", P, iStep, time)
    gmsh.gmshOutput_element(
        gmsh_out, "Euler-Almansi strain", euler, iStep, time)
    gmsh.gmshOutput_element(gmsh_out, "Hencky strain", hencky, iStep, time)
    gmsh.gmshOutput_element(gmsh_out, "Cauchy stress", sigma, iStep, time)

    # break


# print("U[nDofs]")

scale = 100

for i in range(nNodes):
    nodei = mesh.getNode(i)
    nodei.X += scale * U_notflat[i, :dim]

plot.plot_mesh(mesh, label="After", color="r")
plt.legend()
plt.show()


gmsh_out.close()  # close gmsh output file
# ######################
# #%% end of simulation
# #####################
# # print convergence status and solver output (iteration, CPU time)
# print('####################################################')
# if (norm<=test):
#     print('############### SIMULATION COMPLETED ###############')
# else:
#     print('!SIMULATION FAILED: NO CONVERGENCE AT LAST TIME STEP ')
# print('####################################################')
# itermoy /= nSteps
# print('Mean number of iteration per step : %.0f' % itermoy)
# interval = datetime.now() - start_time #measure time for the simulation
# print('Total time of computation:',str(interval))
#
# ####################
# #%% post-processing
# ###################
# print('*-----------------------*')
# print('*--- POST PROCESSING ---*')
# print('*-----------------------*')
# #
# #--- define general plot settings
# #
# import matplotlib as mpl
# mpl.rcParams['font.size'] = 14.0
# mpl.rcParams['lines.linewidth'] = 2
# mpl.rcParams['lines.markersize'] = 8
# mpl.rcParams['lines.color'] = 'r'
# mpl.rcParams['axes.grid'] = True
# #
# #--- plot loading curve (prescribed disp vs. reaction)
# #
# print('Loading curve...')
# # create a new figure and associated axes
# fig1 = plt.figure()
# ax1 = fig1.add_subplot(111)
# # plot internal displacement vs. pressure
# reactions = np.array(reactions) #shape (nSteps x nNodesOnInteriorBoundary)
# displacements = np.array(displacements) #shape (nSteps x nNodesOnInteriorBoundary)
# internal_displacement = np.mean(displacements,axis=1) # all nodal disp. are equal, so the mean = rdisp
# internal_pressure = reactions.sum(axis=1)/(np.pi*Ri/2)
# # compute theoretical result pi = c1*u(Ri) + c2 (cf. MEMCO TD9 or TP3)
# lam = material.prop['1ST_LAME_CONSTANT'] ; mu = material.prop['2ND_LAME_CONSTANT']
# c0 = mu*Ri**3 + (mu+lam)*Re**2*Ri #constant
# c1 = -2*mu*(lam+mu)*(Ri-Re)*(Ri+Re)/c0 #constant
# c2 = mu*Re**2*Ri/c0 #constant
# internal_pressure_th = c1*internal_displacement + c2
# # plot
# ax1.plot(internal_displacement*1e3,internal_pressure_th*1e-6,'-k',label="LinElast theory") # plot in (mm) and (MPa)
# ax1.plot(internal_displacement*1e3,internal_pressure*1e-6,'--o',label="FE") # plot in (mm) and (MPa)
# #configure plot
# ax1.set_xlabel('Prescribed displacement $u_d$ (mm)')
# ax1.set_ylabel('Internal pressure (MPa)')
# ax1.legend() #add a legend
# # savefig to .png file, at the same spot as the .msh file
# fig_out1 = outputfile[:-4]+'_loading-curve.png'
# fig1.tight_layout()
# fig1.savefig(fig_out1,dpi=300,bbox_inches='tight')
# fig1.show()
# #
# #--- plot evolution of the displacement along the bottom line
# #
# print('Displacement...')
# # get nodes number and first coordinates along the bottom line
# bottomLine = np.array([(i,mesh.getNode(i).getX(0)) for i in range(nNodes) if mesh.getNode(i).getX(1)<=Ymin])
# bottomNodes = bottomLine[:,0].astype('int')
# bottomCoordX = bottomLine[:,1]
# # get nodal displacement along the bottom line
# bottomU = U.reshape(nNodes,dim)[bottomNodes][:,0].flatten()
# # compute theoretical result u(r) = a*r + b/r (cf. MEMCO TD9 or TP3)
# lam = material.prop['1ST_LAME_CONSTANT'] ; mu = material.prop['2ND_LAME_CONSTANT']
# pi = internal_pressure[-1] #interior pressure at last time step
# Req = 1/(Re**2-Ri**2) ; A = (Ri**2*pi - Re**2)*Req ; B = Ri**2*Re**2*Req*pi #constants
# a = A/(2*(lam+mu)) ; b = B/(2*mu) # constants
# bottomU_linTH = a * bottomCoordX + b / bottomCoordX # theoretical nodal displacements u_r = a*r + b/r
# # create figure and axis
# fig2 = plt.figure()
# ax2 = fig2.add_subplot(111)
# # sort values according to ascending X-coordinates for the plot
# ind = np.unravel_index(np.argsort(bottomCoordX,axis=None),bottomCoordX.shape)
# # plot
# ax2.plot(bottomCoordX[ind]*1e3,bottomU_linTH[ind]*1e3,'-k',label="LinElast theory")
# ax2.plot(bottomCoordX[ind]*1e3,bottomU[ind]*1e3,'--o',label="FE")
# #configure plot
# ax2.set_xlabel('Coordinate along bottom line $X$ (mm)')
# ax2.set_ylabel('Prescribed displacement $u_d$ (mm)')
# ax2.legend() #add a legend
# # savefig to .png file, at the same spot as the .msh file
# fig_out2 = outputfile[:-4]+'_disp_compare_th.png'
# fig2.tight_layout()
# fig2.savefig(fig_out2,dpi=300,bbox_inches='tight')
# fig2.show()
