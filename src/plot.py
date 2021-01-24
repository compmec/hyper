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


def plot_mesh(mesh, label=None, color=None):

    if label is None:
        label = "Mesh"
    if color is None:
        color = "k"

    show_nodename = False

    nElements = mesh.getNElements()
    nNodes = mesh.getNNodes()
    X_all = np.zeros(nNodes)
    Y_all = np.zeros(nNodes)
    all_connections = []
    for i in range(nElements):
        elem = mesh.getElement(i)
        connections = elem.getNodes()
        all_connections.append(connections)
    if max([len(c) for c in all_connections]) == 4:
        return
    for i in range(nNodes):
        Nodei = mesh.getNode(i)
        X_all[i] = Nodei.getX(0)
        Y_all[i] = Nodei.getX(1)

    # print("X_all = ")
    # print(X_all)
    # print("Y_all = ")
    # print(Y_all)
    # print("all_connections = ")
    # print(all_connections)

    if show_nodename:
        for i in range(nNodes):
            txt = str(mesh.getNode(i).id)
            txt = str(i)
            x = X_all[i]
            y = Y_all[i]
            plt.annotate(txt, (x, y))

    plt.triplot(X_all, Y_all, all_connections, color=color, label=label)
    plt.axis("equal")
    # plt.show()


if __name__ == "__main__":
    # -- read mesh from Gmsh file
    #
    msh_folder = "../msh/"
    # basename = "triangle.msh"
    basename = "cylinder-custom.msh"
    inputfnam = msh_folder + basename
    # inputfnam = msh_folder + "meshfile-tri8.msh"
    meshfile = open(inputfnam, 'r')
    mesh = gmsh.gmshInput_mesh(meshfile)
    meshfile.close()

    nNodes = mesh.getNNodes()
    nEleme = mesh.getNElements()

    print("Read mesh with %d nodes and %d elements" % (nNodes, nEleme))

    plot_mesh(mesh)
    plt.show()
