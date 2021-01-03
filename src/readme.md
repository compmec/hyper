# Source

In this folder there are the main files used to make calculus and solve the problem.


### Files

* The file ```tensor.py``` contains all the tensor operations that we need to solve the problem. It uses ```numpy``` as base library.
* The file ```gmsh2.py``` are responsible for read and writing in the ```.msh```files. As there is a ```gmsh``` python's native library, it could be a problem, so we changed the name. 
* The file ```mesh.py``` contains the classes to store the informations of the mesh
    + ```Node``` contains the position of the points
    + ```Element``` contain the relation of many ```Node```
    + ```MeshData``` storage all the nodes and the elements
* The file ```geometry.py``` that contains all the functions of interpolations and the integrations points. 
    + The functions of interpolation is explained in [Base Element's Geometry](https://github.com/carlos-adir/Non-linear-mechanics/wiki/Base-Elements'-Geometry)
        - Elements T3, T6, Q4 and Q8
    + The integrations points are shown in [Gauss Integration Points](https://github.com/carlos-adir/Non-linear-mechanics/wiki/Gauss-Integration-Points)
        - In a triangle and in a square.
        - Order ```O(1)```, ```O(h)``` and ```O(h^2)```
* The file ```fem.py``` has the elements to calculate the values of the tensors ```F``` and ```E``` if we already know the value of the displacement field ```u``` of the elements
    + ```FiniteElement``` class stores and treat the informations of only one element
    + ```FEModel``` class stores and treat the information of all elements
* The file ```elasticity.py``` contains the hyperelastic models of the material.
    + ```StVenantKirchhoffElasticity``` class is the [Saint Venant-Kirchhoff][stvenantkirchhoff] model
    + ```NeoHookeanElasticity``` class is the [NeoHookean][neohookean] model
* The file ```plot.py``` contains the functions to plot the mesh. It's an alternative to visualize the data instead using the software [Gmsh][gmsh_website]. It's not yet finished.

### Architecture

The relation of import and dependency between the files are:

<img src="https://raw.githubusercontent.com/carlos-adir/Non-linear-mechanics/docs/img/architecturesource.png" alt="architecture" width="600"/>


[stvenantkirchhoff]: https://en.wikipedia.org/w/index.php?title=Hyperelastic_material&oldid=993354665
[neohookean]: https://en.wikipedia.org/w/index.php?title=Neo-Hookean_solid&oldid=980304435
[gmsh_website]: https://gmsh.info/