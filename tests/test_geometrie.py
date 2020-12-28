# -*- coding: utf-8 -*-


try:
    from test_unit import TU
except ModuleNotFoundError:
    print("Main file to test was not found: test_unit.py")

try:
    TU.include_path("src")
    import geometrie
except ModuleNotFoundError:
    print("Import file not found: tensor.py")


class TUGeometrie(TU):
    """
    Class for testing the geometrie module
    """

    # FUNCTIONS = ["T3", "T6", "Q4", "Q8"]
    FUNCTIONS = ["T3"]

    def __init__(self):
        super(TUGeometrie, self).__init__()

    @staticmethod
    def test_interpolation(points, N):
        n = len(points)
        for i in range(n):
            Np = N(points[i])
            for j in range(n):
                if i == j:
                    if Np[j] != 1:
                        raise Exception("N(pi)[i] != 1 ---> Impossible")
                else:
                    if Np[j] != 0:
                        raise Exception(
                            "N(pi)[j] != 0 ---> Impossible")

    @staticmethod
    def test_derivative(points, dN):
        n = len(points)
        for i in range(n):
            dNp = dN(points[i])
            for j in range(n):
                pass
                # if i == j:
                #     if dNp[i] != 1:
                #         raise Exception("dN(pi)[i] != 1 ---> Impossible")
                # else:
                #     if dNp[i] != 0:
                #         raise Exception("dN(pi)[j] != 0 ---> Impossible")

    @staticmethod
    def T3():
        N = geometrie.SFT3.N
        dN = geometrie.SFT3.dN
        P = geometrie.SFT3.P
        points = P()
        interp = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
        deriva = [((-1, -1), (1, 0), (0, 1)),
                  ((-1, -1), (1, 0), (0, 1)),
                  ((-1, -1), (1, 0), (0, 1))]
        TUGeometrie.test_interpolation(points, N)
        TUGeometrie.test_derivative(points, dN)

        # n = len(points)
        # for i, p in enumerate(points):
        #     Np = N(p)
        #     for j in range(n):
        #         if Np[j] != interp[i][j]:
        #             raise Exception("Not interpolation")
        # for i, p in enumerate(points):
        #     dNp = dN(p)
        #     for j in range(n):
        #         for k in range(2):
        #             if dNp[j][k] != deriva[i][j][k]:
        #                 raise Exception("Not derivation")


if __name__ == "__main__":
    test = TUGeometrie()
    test.run()
