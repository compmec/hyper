[![Build Status][build-img]][build-url]

# Hyper-elastic materials

### Introduction

This repository contains the code used in the cours HYPER (Numerical Methods for NonLinear Mechanics) at the [École Centrale de Nantes](https://www.ec-nantes.fr/).

This cours describes the fundamentals of Hyper-elastic materials and how to implement finite element method with them.


### Subject

Hyper-elastic materials are used when the hypothesis of small deformations is not valid anymore.

For exemple, applying pressure in a steel cylinder is different from applying pressure in a rubber pipe, which makes the pipe inflates and reduces its thickness:

![Neo Hook model - Rubber](https://raw.githubusercontent.com/carlos-adir/Non-linear-mechanics/docs/img/LinearMechanics.gif)

The main difference between the approachs are:

* Linear elastic materials
    + Valid for small deformations
    + The stiffness matrix is constant
    + Linear problem: solves a directly a single linear system
* Hyper-elastic materials
    + Valid for small and large deformations
    + Needs recalculation of stiffness matrix, to consider the deformed space
    + Non-linear problems: Needs to 'walk', solving a linear system for each step

![Difference linear / nonlinear](https://ars.els-cdn.com/content/image/3-s2.0-B9780323898225000104-f04-01-9780323898225.jpg)


### How to use

Python is used to implement the finite element method, along with the libraries:

* [Numpy](https://numpy.org/doc/): Used for tensor calculs
* [Gmsh](https://gmsh.info/): Used to get the mesh and elements
* [Scipy](https://scipy.org/): Used to get the mesh and elements

To install them, you can either use [Anaconda](https://www.anaconda.com/) that installs all you need, or you can install them manually:

```
pip install numpy
pip install scipy
pip install gmsh
```

After getting the libraries and downloading your own copy, you can test it with `pytest`:

```
pip install pytest
pytest
```

There are some examples in the folder `examples` explaining how to run


<!-- Badges: -->

[build-img]: https://github.com/compmec/hyper/actions/workflows/build.yaml/badge.svg
[build-url]: https://github.com/compmec/hyper/actions/workflows/build.yaml
