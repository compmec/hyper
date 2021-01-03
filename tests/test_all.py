from test_unit import TU
from test_tensor import TUTensor
from test_mesh import TUMesh
from test_geometry import TUGeometry
from test_fem import TUFem
from test_elasticity import TUElasticity
from test_plot import TUPlot
from test_residual import TUResidual


if __name__ == "__main__":
    TU.VERBOSE = True
    test = TUTensor()
    test.run()
    test = TUMesh()
    test.run()
    test = TUGeometry()
    test.run()
    test = TUPlot()
    test.run()
    test = TUFem()
    test.run()
    test = TUElasticity()
    test.run()
    test = TUResidual()
    test.run()
