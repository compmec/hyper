from test_unit import TU
from test_tensor import TUTensor
from test_mesh import TUMesh
from test_geometry import TUGeometry
from test_fem import TUFem
from test_elasticity import TUElasticity
from test_plot import TUPlot
from test_residual import TUResidual


if __name__ == "__main__":
    TU.VERBOSE = False
    TUClasses = [TUTensor, TUMesh, TUGeometry,
                 TUPlot, TUFem, TUElasticity, TUResidual, test_tangent]
    counter = 0
    for TUclass in TUClasses:
        test = TUclass()
        result = test.run()
        if result == TU.SUCESS:
            counter += 1
    print("Sucess: [%d/%d] " % (counter, len(TUClasses)))
