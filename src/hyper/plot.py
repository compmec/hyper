import numpy as np
from matplotlib import pyplot as plt


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

    if show_nodename:
        for i in range(nNodes):
            txt = str(mesh.getNode(i).id)
            txt = str(i)
            x = X_all[i]
            y = Y_all[i]
            plt.annotate(txt, (x, y))

    plt.triplot(X_all, Y_all, all_connections, color=color, label=label)
    plt.axis("equal")

