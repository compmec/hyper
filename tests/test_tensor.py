# -*- coding: utf-8 -*-


try:
    import numpy as np
    from numpy import linalg as la
except ModuleNotFoundError:
    print("The numpy library was not found")


try:
    from test_unit import TU
except ModuleNotFoundError:
    print("Main file to test was not found: test_unit.py")

try:
    TU.include_path("src")
    import tensor
except ModuleNotFoundError:
    print("Import file not found: tensor.py")


class TUTensor(TU):
    """
    Class for testing the tensor module
    """

    FUNCTIONS = ["vector",
                 "tensor",
                 "I",
                 "outerProd",
                 "det",
                 "trace",
                 "inv",
                 "rightCauchyGreen",
                 "leftCauchyGreen",
                 "tensor4",
                 "II",
                 "IISym",
                 "KK",
                 "outerProd4"]

    def __init__(self):
        super(TUTensor, self).__init__()

    @staticmethod
    def vector():
        for number in range(3):
            if number == 0:
                tensor.vector()
            elif number == 1:
                my_vector = tensor.vector()
                my_vector[0] = 1.0
                my_vector[1] = 2.0
            elif number == 2:
                my_vector = tensor.vector(3)
                my_vector[0] = 1.5
                my_vector[2] = -3.0

    @staticmethod
    def tensor():
        tensor.tensor()

    @staticmethod
    def I():
        tensor.I()

    @staticmethod
    def outerProd():
        v = np.array([1, 2])
        output = tensor.outerProd(v, v)

        expected = np.array([[1, 2],
                             [2, 4]])
        if la.norm(output - expected) > 1e-14:
            raise Exception("Not equal")

    @staticmethod
    def det():
        for number in range(3):
            if number == 0:
                B = np.array([[1, 2],
                              [2, 4]])
                output = tensor.det(B)

                expected = 0
            elif number == 1:
                B = np.array([[1, 2],
                              [2, 5]])
                output = tensor.det(B)

                expected = 1
            if la.norm(output - expected) > 1e-14:
                raise Exception("Not equal")

    @staticmethod
    def trace():
        B = np.array([[1, 2],
                      [2, 4]])
        output = tensor.trace(B)
        expected = 5
        if la.norm(output - expected) > 1e-14:
            raise Exception("Not equal")

    @staticmethod
    def inv():
        for number in range(2):
            if number == 0:
                F = np.array([[2, 3],
                              [2, 5]])

                expected = np.array([[1.25, -0.75],
                                     [-0.5, 0.5]])
            elif number == 1:
                F = np.array([[1, 0],
                              [0, 1]])
                expected = np.array([[1, 0],
                                     [0, 1]])
            output = tensor.inv(F)
            if la.norm(output - expected) > 1e-14:
                raise Exception("Not equal")

    @staticmethod
    def rightCauchyGreen():
        F = np.array([[2, 3],
                      [2, 5]])
        output = tensor.rightCauchyGreen(F)
        expected = np.array([[8, 16],
                             [16, 34]])
        if la.norm(output - expected) > 1e-14:
            raise Exception("Not equal")

    @staticmethod
    def leftCauchyGreen():
        F = np.array([[2, 3],
                      [2, 5]])
        output = tensor.leftCauchyGreen(F)
        expected = np.array([[13, 19],
                             [19, 29]])
        if la.norm(output - expected) > 1e-14:
            raise Exception("Not equal")

    @staticmethod
    def tensor4():
        output = tensor.tensor4()
        expected = np.zeros((2, 2, 2, 2))
        # expected[0, 0, 0, 0] = 1
        if la.norm(output - expected) > 1e-14:
            raise Exception("Not equal")

    @staticmethod
    def II():
        output = tensor.II()
        expected = np.array([[[[1, 0], [0, 0]], [[0, 1], [0, 0]]], [
                            [[0, 0], [1, 0]], [[0, 0], [0, 1]]]])
        if la.norm(output - expected) > 1e-14:
            raise Exception("Not equal")

    @staticmethod
    def IISym():
        output = tensor.IISym()
        expected = np.array([[[[1,   0],  [0,   0]],  [[0,   0.5],  [0.5, 0]]],               [
                            [[0,   0.5],  [0.5, 0]],  [[0,   0],  [0,   1.]]]])
        if la.norm(output - expected) > 1e-14:
            raise Exception("Not equal")

    @staticmethod
    def KK():
        tensor.KK()

    @staticmethod
    def outerProd4():
        C = np.array([[8, 16],
                      [16, 34]])
        output = tensor.outerProd4(C, C)
        expected = [[[[64, 128], [128, 272]], [[128, 256], [256, 544]]], [
            [[128, 256], [256, 544]], [[272, 544], [544, 1156]]]]
        if la.norm(output - expected) > 1e-14:
            raise Exception("Not equal")


if __name__ == "__main__":
    test = TUTensor()
    test.run()
