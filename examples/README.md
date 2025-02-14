# Example

This folder contains examples of how to use the code made.


### Cylinder

For the moment, there's only one file: **cylinder.py**

Just run it on your terminal:

```
python cylinder.py
```

The result will depend on the used parameters, which are:

* material: meaning the stiffness and the used law
    + **StVenant: steel**
    + NeoHookean: rubber
* mesh: the definition of nodes position and elements (triangles, squares):
    + **cylinder-custom.msh**
    + cylinder-tri_coarse.msh
    + cylinder-tri_refined.msh
    + cylinder-quad_coarse.msh
    + cylinder-quad_refined.msh
* iteration: tells how fast we want to approach the solution 
    + nSteps: more steps means solving more linear systems, but increases change of convergence
    + rdisp: the charge applied on the center
* solver: tells how we want to solve one step, from a time (t0) to a time (t0 + delta T)
    + precision: the zero machine
    + tolerance: tell when should we stop solving the linear system
    + itermax: the number maximal of iteration to solve the linear system. The computation is stoped if there's no convergence


The default is marked in **bold** and you can change it as you want:

```
============================================
=                                          =
=      2D plane strain thick cylinder      =
=                                          =
============================================
############################
#         MATERIAL         #
############################
    Material: steel
           E = 2.10e+07
           v = 3.00e-01
       model = StVenant

############################
#   ITERATION PARAMETERS   #
############################
     Number of steps = 10
             delta T = 0.1
 radial displacement = 1.00e-04

############################
#     SOLVER PARAMETERS    #
############################
    PRECISION = 1.00e-15
    TOLERANCE = 1.00e-06
      ITERMAX = 10.00

############################
#     FILES AND FOLDERS    #
############################
    mesh = cylinder-custom.msh
    geom = cylinder.geo
    Geometry cylinder:
        Internal radius Ri = 2.00e-02 m
        External radius Re = 6.00e-02 m
    Mesh read:
        File: msh/cylinder-custom.msh
        Nodes: 8
        Elements: 8

############################
#      PRE PROCESSING      #
############################
    Mesh information:
        Total nodes: 8
        Total elements: 8
        Dimention: 2
    Number of nodes:
        Left: 2
        Center: 3
        Bottom: 2
    Degrees of freedom:
        Total: 16
        BC Left: 2
        BC Center: 6
        BC Bottom: 2
        BC Total: 8
        Free: 8

############################
#        SIMULATION        #
############################
--> Step [1/10]
    time = 0.10/1.00
    disp = 1.00e-05
    *** Iteration 00: residual norm =2.465e+02 ***
    *** Iteration 01: residual norm =6.998e-02 ***
    *** Iteration 02: residual norm =7.345e-09 ***
    COMPLETE: norm < test
        norm = 7.34500e-09
        test = 2.46498e-04
...
--> Step [10/10]
    time = 1.00/1.00
    disp = 1.00e-04
    *** Iteration 00: residual norm =2.572e+02 ***
    *** Iteration 01: residual norm =5.599e+00 ***
    *** Iteration 02: residual norm =4.787e-05 ***
    COMPLETE: norm < test
        norm = 4.78701e-05
        test = 2.57201e-04
```


**PS:** The geo folder contains the initial and vectorized geometry. It's used to make the `.msh` files, but in the simulation it has no much importance.
It's used only to get the value of the internal radius, to then get the points near the inner center to apply the displacement.
