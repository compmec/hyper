# Source files

In this folder there are the files used to make calculs and the tests files.


### Library files

* The file ```gmsh.py``` are responsable for read and writing in the ```.msh```files.
* The file ```mesh.py``` contains the classes to store the informations of the mesh
    + ```Node``` contains the position of the points
    + ```Element``` contain the relation of many ```Node```
    + ```MeshData``` storage all the nodes and the elements
* The file ```geometrie.py``` that contains all the functions of interpolations and the integrations points. 
    + The functions of interpolation is explained in [Base Element's Geometry](https://github.com/carlos-adir/Non-linear-mechanics/wiki/Base-Elements'-Geometry)
    + The integrations points are shown in [Gauss Integration Points](https://github.com/carlos-adir/Non-linear-mechanics/wiki/Gauss-Integration-Points)
* The file ```fem.py``` has the elements to calculate the values of the tensors ```F``` and ```E``` if we already know the value of the displacement field ```u``` of the elements
    + ```FiniteElement``` stores and treat the informations of only one element
    + ```FEModel``` stores and treat the information of all elements


### Tests files

* The file ```test_tensor.py``` tests the functions in the file ```tensor.py```
* The file ```test_mesh.py``` tests the implementation of the file ```mesh.py```, we are assumming here that we don't need to test the file ```gmsh.py```
* The file ```test_fem.py``` tests the file ```fem.py``` and the file ```geometrie.py``` 