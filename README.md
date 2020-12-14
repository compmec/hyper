# Non-linear-mechanics

### Introduction

This repository contains the code used in the cours HYPER(Numerical Methods for NonLinear Mechanics) at the [École Centrale de Nantes](https://www.ec-nantes.fr/).


### Subject

The non-linear mechanics is used when the hypothesis of small deformations is not valid. For exemple, we have the images below that shows when it's applied a big force in the center.


![](https://raw.githubusercontent.com/carlos-adir/Non-linear-mechanics/docs/img/LinearMechanics.png)

The main difference between the approachs are:

* Linear mechanics
    + Integration doesn't consider the deformed space
    + Eulerien's description
    + Cauchy's strain tensor
* Non-linear mechanics
    + Integration considers the deformed space
    + Lagrangian's description
    + Green-Lagrangian's strain tensor



### Coding and librarys

For implementation, we use Python with the libraries:

* [Numpy](https://numpy.org/doc/): Used for tensor calculs
* [Gmsh](https://gmsh.info/): Used to get the mesh and elements

To use the codes, you just need these packages and the python installed. The easiest way to do it is using the [Anaconda](https://www.anaconda.com/) that installs everything you need.


### Documentation

All documentation is kept in [our wiki](https://github.com/carlos-adir/Non-linear-mechanics/wiki), with usage and examples.


### Authors

* Carlos Adir Ely Murussi Leite
* Clément Fuchs