try:
    import numpy as np
    from matplotlib import pyplot as plt
except ModuleNotFoundError:
    print("The numpy library was not found")


try:
    import gmsh2 as gmsh
    import geometry
except ModuleNotFoundError:
    print("Import file not found: gmsh2.py")


def plot_mesh(mesh):
    nElements = mesh.getNElements()
    nNodes = mesh.getNNodes()
    X_all = np.zeros(nNodes)
    Y_all = np.zeros(nNodes)
    all_connections = []
    for i in range(nElements):
        elem = mesh.getElement(i)
        connections = elem.getNodes()
        all_connections.append(connections)
    for i in range(nNodes):
        Nodei = mesh.getNode(i)
        X_all[i] = Nodei.getX(0)
        Y_all[i] = Nodei.getX(1)

    plt.triplot(X_all, Y_all, all_connections, color="k")
    plt.axis("equal")
    plt.title("Mesh")
    plt.show()


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

    plot_mesh(mesh)
