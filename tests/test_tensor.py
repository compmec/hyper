# -*- coding: utf-8 -*-
import numpy as np

from hyper import tensor


def test_vector():
    # Testing creation
    vectest = tensor.vector()
    assert type(vectest) is np.ndarray

    # Testing attribution of individual values
    vectest[0] = 1.0
    vectest[1] = 2.0
    assert type(vectest) is np.ndarray
    assert vectest.ndim == 1
    assert vectest.shape == (2,)
    np.testing.assert_array_equal(vectest, [1, 2])

    # Testing attribution of array values
    vectest[:] = [3, 4]
    assert type(vectest) is np.ndarray
    assert vectest.ndim == 1
    assert vectest.shape == (2,)
    np.testing.assert_array_equal(vectest, [3, 4])

    # Testing dim = 3 array
    vectest = tensor.vector(3)
    vectest[:] = [-1, 2, 1]
    assert type(vectest) is np.ndarray
    assert vectest.ndim == 1
    assert vectest.shape == (3,)
    np.testing.assert_array_equal(vectest, [-1, 2, 1])


def test_tensor():
    # Testing tensor with dim = 2
    Mtest = tensor.tensor()
    assert type(Mtest) is np.ndarray
    assert Mtest.ndim == 2
    assert Mtest.shape == (2, 2)

    # Testing tensor with dim = 3
    Mtest = tensor.tensor(3)
    assert type(Mtest) is np.ndarray
    assert Mtest.ndim == 2
    np.testing.assert_array_equal(Mtest.shape, (3, 3))


def test_I():
    Itest = tensor.I()
    Igood = [[1, 0], [0, 1]]
    assert type(Itest) is np.ndarray
    assert Itest.shape == (2, 2)
    np.testing.assert_array_equal(Itest, Igood)

    Itest = tensor.I(3)
    Igood = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    assert type(Itest) is np.ndarray
    assert Itest.shape == (3, 3)
    np.testing.assert_array_equal(Itest, Igood)


def test_outerProd():
    v = [1, 2]
    Mgood = [[1, 2], [2, 4]]
    Mtest = tensor.outerProd(v, v)
    assert Mtest.ndim == 2
    assert Mtest.shape == (2, 2)
    np.testing.assert_array_equal(Mtest, Mgood)


def test_det():
    M = [[1, 2], [2, 4]]
    detgood = 0
    dettest = tensor.det(M)
    assert dettest == detgood

    M = [[1, 2], [2, 5]]
    detgood = 1
    dettest = tensor.det(M)
    assert dettest == detgood


def test_trace():
    M = [[1, 2], [2, 4]]
    trgood = 5
    trtest = tensor.trace(M)
    assert trtest == trgood

    M = [[1, 2], [2, 5]]
    trgood = 6
    trtest = tensor.trace(M)
    assert trtest == trgood


def test_inv():
    M = [[1, 0], [0, 1]]
    Vgood = [[1, 0], [0, 1]]
    Vtest = tensor.inv(M)
    assert Vtest.ndim == 2
    assert Vtest.shape == (2, 2)
    np.testing.assert_array_equal(Vtest, Vgood)

    M = [[2, 3], [2, 5]]
    Vgood = [[1.25, -0.75], [-0.5, 0.5]]
    Vtest = tensor.inv(M)
    assert Vtest.ndim == 2
    assert Vtest.shape == (2, 2)
    np.testing.assert_array_equal(Vtest, Vgood)


def test_rightCauchyGreen():
    F = [[2, 3], [2, 5]]
    Cgood = [[8, 16], [16, 34]]
    Ctest = tensor.rightCauchyGreen(F)
    assert Ctest.ndim == 2
    assert Ctest.shape == (2, 2)
    np.testing.assert_array_equal(Ctest, Cgood)


def test_leftCauchyGreen():
    F = [[2, 3], [2, 5]]
    Bgood = [[13, 19], [19, 29]]
    Btest = tensor.leftCauchyGreen(F)
    assert Btest.ndim == 2
    assert Btest.shape == (2, 2)
    np.testing.assert_array_equal(Btest, Bgood)


def test_tensor4():
    Mtest = tensor.tensor4()
    assert Mtest.ndim == 4
    assert Mtest.shape == (2, 2, 2, 2)
    assert np.sum(Mtest) == 0

    Mtest = tensor.tensor4(3)
    assert Mtest.ndim == 4
    assert Mtest.shape == (3, 3, 3, 3)
    assert np.sum(Mtest) == 0


def test_II():
    IIgood = [
        [[[1, 0], [0, 0]], [[0, 1], [0, 0]]],
        [[[0, 0], [1, 0]], [[0, 0], [0, 1]]],
    ]
    IItest = tensor.II()
    assert IItest.ndim == 4
    assert IItest.shape == (2, 2, 2, 2)
    np.testing.assert_array_equal(IItest, IIgood)


def test_IISym():
    IItest = tensor.IISym()
    IIgood = [
        [[[1, 0], [0, 0]], [[0, 0.5], [0.5, 0]]],
        [[[0, 0.5], [0.5, 0]], [[0, 0], [0, 1]]],
    ]
    assert IItest.ndim == 4
    assert IItest.shape == (2, 2, 2, 2)
    np.testing.assert_array_equal(IItest, IIgood)


def test_KK():
    KKtest = tensor.KK()
    KKgood = [
        [[[1, 0], [0, 1]], [[0, 0], [0, 0]]],
        [[[0, 0], [0, 0]], [[1, 0], [0, 1]]],
    ]
    assert KKtest.ndim == 4
    assert KKtest.shape == (2, 2, 2, 2)
    np.testing.assert_array_equal(KKtest, KKgood)


def test_outerProd4():
    C = [[8, 16], [16, 34]]
    Mgood = [
        [[[64, 128], [128, 272]], [[128, 256], [256, 544]]],
        [[[128, 256], [256, 544]], [[272, 544], [544, 1156]]],
    ]
    Mtest = tensor.outerProd4(C, C)
    assert Mtest.ndim == 4
    assert Mtest.shape == (2, 2, 2, 2)
    np.testing.assert_array_equal(Mtest, Mgood)
