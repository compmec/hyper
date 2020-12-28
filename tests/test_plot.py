# -*- coding: utf-8 -*-
#


try:
    from test_unit import TU
except ModuleNotFoundError:
    print("Main file to test was not found: test_unit.py")

try:
    TU.include_path("src")
    import gmsh2 as gmsh
    import plot
except ModuleNotFoundError:
    print("Import file not found: tensor.py")


class TUPlot(TU):

    def __init__(self):
        super(TUPlot, self).__init__()


if __name__ == "__main__":
    # -- read mesh from Gmsh file
    #
    msh_folder = "../msh/"
    inputfnam = msh_folder + "triangle-tri56.msh"
    # inputfnam = msh_folder + "triangle-quad44.msh"
    meshfile = open(inputfnam, 'r')
    mesh = gmsh.gmshInput_mesh(meshfile)
    meshfile.close()

    nNodes = mesh.getNNodes()
    nEleme = mesh.getNElements()
    print("Read mesh with %d nodes and %d elements" % (nNodes, nEleme))

    plot.plot_mesh(mesh)
