# Tests

In this folder there are the files used test the main code.
With them, it's easier to test the code and debug if something happens.


### Files

There are two different files that we have to give attention:

* ```test_unit.py``` contains the basic class ```TU```(Test Unit) that the others test classes inherit from.
* ```test_all.py``` contains all the tests, it's used if we want to test all the tests at the same time, instead running each test file separately.

The others files are to test a specific part/module:

1. The file ```test_tensor.py``` tests the file ```src/tensor.py```.
2. The file ```test_mesh.py``` tests the file ```src/mesh.py```
3. The file ```test_geometry.py``` test the file ```src/geometry.py```
4. The file ```test_plot.py``` test the file ```src/plot.py``` 
5. The file ```test_fem.py``` tests the file ```src/fem.py``` that uses the ```src/geometry.py```. For example, create the finite element and integrate in the geometry.
6. The files ```test_elasticity.py``` tests the computation of piola's tensor using the hyperelastic materials in the ```src/elasticity.py```
7. The file ```test_residual.py``` tests the part of calculating the internal forces in the file ```src/fem.py```.



[stvenantkirchhoff]: https://en.wikipedia.org/w/index.php?title=Hyperelastic_material&oldid=993354665
[neohookean]: https://en.wikipedia.org/w/index.php?title=Neo-Hookean_solid&oldid=980304435
[gmsh_website]: https://gmsh.info/
