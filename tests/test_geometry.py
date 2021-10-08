# -*- coding: utf-8 -*-

import pytest
import numpy as np
from hyper import geometry


def getRandomPointTriangle():
    # We want points (x, y) that satisfy
    # 0 <= x <= 1, 0 <= y <= 1 and (x+y) <= 1
    x = np.random.rand(1)
    y = np.random.rand(1) * (1 - x)
    return (x, y)


def getRandomPointSquare():
    # We want points (x, y) that satisfy
    # -1 <= x <= 1, -1 <= y <= 1
    x = 2 * np.random.rand(1) - 1
    y = 2 * np.random.rand(1) - 1
    return (x, y)


def test_SumValuesSFT3():
    Ntests = 100
    for i in range(Ntests):
        point = getRandomPointTriangle()
        Np = geometry.SFT3.N(point)
        assert np.sum(Np) == pytest.approx(1)


def test_SumValuesSFT6():
    Ntests = 100
    for i in range(Ntests):
        point = getRandomPointTriangle()
        Np = geometry.SFT6.N(point)
        assert np.sum(Np) == pytest.approx(1)


def test_SumValuesSFQ4():
    Ntests = 100
    for i in range(Ntests):
        point = getRandomPointTriangle()
        Np = geometry.SFQ4.N(point)
        assert np.sum(Np) == pytest.approx(1)


def test_SumValuesSFQ8():
    Ntests = 100
    for i in range(Ntests):
        point = getRandomPointTriangle()
        Np = geometry.SFQ8.N(point)
        assert np.sum(Np) == pytest.approx(1)


def test_InterpolationSFT3():
    Ps = geometry.SFT3.P()  # Element's point
    n = len(Ps)  # Number of points
    for i, Pi in enumerate(Ps):
        V = geometry.SFT3.N(Pi)
        for j in range(n):
            if i == j:
                assert V[j] == 1
            else:
                assert V[j] == 0


def test_InterpolationSFT6():
    Ps = geometry.SFT6.P()  # Element's point
    n = len(Ps)  # Number of points
    for i, Pi in enumerate(Ps):
        V = geometry.SFT6.N(Pi)
        for j in range(n):
            if i == j:
                assert V[j] == 1
            else:
                assert V[j] == 0


def test_InterpolationSFQ4():
    Ps = geometry.SFQ4.P()  # Element's point
    n = len(Ps)  # Number of points
    for i, Pi in enumerate(Ps):
        V = geometry.SFQ4.N(Pi)
        for j in range(n):
            if i == j:
                assert V[j] == 1
            else:
                assert V[j] == 0


def test_InterpolationSFQ8():
    Ps = geometry.SFQ8.P()  # Element's point
    n = len(Ps)  # Number of points
    for i, Pi in enumerate(Ps):
        V = geometry.SFQ8.N(Pi)
        for j in range(n):
            if i == j:
                assert V[j] == 1
            else:
                assert V[j] == 0


def test_derivativeSFT3():
    # For this function we have constant fields
    # So, doesn't matter the points, it's always the same
    # We test if it changes depending of the point
    Ntests = 100
    dNgood = [[-1, -1], [1, 0], [0, 1]]
    for i in range(Ntests):
        point = getRandomPointTriangle()
        dNtest = geometry.SFT3.dN(point)
        assert dNtest == dNgood


def test_derivativeSFT6():
    P = (0, 0)
    dNgood = [[-3, -3], [4, 0], [-1, 0], [0, 0], [0, -1], [0, 4]]
    dNtest = geometry.SFT6.dN(P)
    assert dNtest == dNgood

    P = (0.5, 0)
    dNgood = [[-1.0, -1.0], [0.0, -2.0], [1.0, 0], [0, 2.0], [0, -1], [0, 2.0]]
    dNtest = geometry.SFT6.dN(P)
    assert dNtest == dNgood

    P = (1, 0)
    dNgood = [[1, 1], [-4, -4], [3, 0], [0, 4], [0, -1], [0, 0]]
    dNtest = geometry.SFT6.dN(P)
    assert dNtest == dNgood

    P = (0, 0.5)
    dNgood = [[-1.0, -1.0], [2.0, 0], [-1, 0], [2.0, 0], [0, 1.0], [-2.0, 0.0]]
    dNtest = geometry.SFT6.dN(P)
    assert dNtest == dNgood

    P = (0.5, 0.5)
    dNgood = [[1.0, 1.0], [-2.0, -2.0], [1.0, 0],
              [2.0, 2.0], [0, 1.0], [-2.0, -2.0]]
    dNtest = geometry.SFT6.dN(P)
    assert dNtest == dNgood

    P = (0, 1)
    dNgood = [[1, 1], [0, 0], [-1, 0], [4, 0], [0, 3], [-4, -4]]
    dNtest = geometry.SFT6.dN(P)
    assert dNtest == dNgood

    P = (1 / 3, 1 / 3)
    dNgood = [[-1 / 3, -1 / 3], [0, -4 / 3],
              [1 / 3, 0], [4 / 3, 4 / 3],
              [0, 1 / 3], [-4 / 3, 0]]
    dNtest = geometry.SFT6.dN(P)
    np.testing.assert_almost_equal(dNtest, dNgood)


def test_derivativeSFQ4():
    P = (-1, -1)
    dNgood = [[-0.5, -0.5], [0.5, 0], [0, 0], [0, 0.5]]
    dNtest = geometry.SFQ4.dN(P)
    assert dNtest == dNgood

    P = (0, -1)
    dNgood = [[-0.5, -0.25], [0.5, -0.25], [0, 0.25], [0, 0.25]]
    dNtest = geometry.SFQ4.dN(P)
    assert dNtest == dNgood

    P = (1, -1)
    dNgood = [[-0.5, 0], [0.5, -0.5], [0, 0.5], [0, 0]]
    dNtest = geometry.SFQ4.dN(P)
    assert dNtest == dNgood

    P = (-1, 0)
    dNgood = [[-0.25, -0.5], [0.25, 0], [0.25, 0], [-0.25, 0.5]]
    dNtest = geometry.SFQ4.dN(P)
    assert dNtest == dNgood

    P = (0, 0)
    dNgood = [[-0.25, -0.25], [0.25, -0.25], [0.25, 0.25], [-0.25, 0.25]]
    dNtest = geometry.SFQ4.dN(P)
    assert dNtest == dNgood

    P = (1, 0)
    dNgood = [[-0.25, 0], [0.25, -0.5], [0.25, 0.5], [-0.25, 0]]
    dNtest = geometry.SFQ4.dN(P)
    assert dNtest == dNgood

    P = (-1, 1)
    dNgood = [[0, -0.5], [0, 0], [0.5, 0], [-0.5, 0.5]]
    dNtest = geometry.SFQ4.dN(P)
    assert dNtest == dNgood

    P = (0, 1)
    dNgood = [[0, -0.25], [0, -0.25], [0.5, 0.25], [-0.5, 0.25]]
    dNtest = geometry.SFQ4.dN(P)
    assert dNtest == dNgood

    P = (1, 1)
    dNgood = [[0, 0], [0, -0.5], [0.5, 0.5], [-0.5, 0]]
    dNtest = geometry.SFQ4.dN(P)
    assert dNtest == dNgood


def test_derivativeSFQ8():
    P = (-1, -1)
    dNtest = geometry.SFQ8.dN(P)
    dNgood = [[-1.5, -1.5], [2, 0],
              [-0.5, 0], [0, 0],
              [0, 0], [0, 0],
              [0, -0.5], [0, 2]]
    assert dNtest == dNgood

    P = (0, -1)
    dNtest = geometry.SFQ8.dN(P)
    dNgood = [[-0.5, -0.5], [0, -0.5],
              [0.5, -0.5], [0, 1],
              [0, -0.5], [0, 0.5],
              [0, -0.5], [0, 1]]
    assert dNtest == dNgood

    P = (1, -1)
    dNtest = geometry.SFQ8.dN(P)
    dNgood = [[0.5, 0], [-2, 0],
              [1.5, -1.5], [0, 2],
              [0, -0.5], [0, 0],
              [0, 0], [0, 0]]
    assert dNtest == dNgood

    P = (-1, 0)
    dNtest = geometry.SFQ8.dN(P)
    dNgood = [[-0.5, -0.5], [1, 0],
              [-0.5, 0], [0.5, 0],
              [-0.5, 0], [1, 0],
              [-0.5, 0.5], [-0.5, 0]]
    assert dNtest == dNgood

    P = (0, 0)
    dNtest = geometry.SFQ8.dN(P)
    dNgood = [[0, 0], [0, -0.5],
              [0, 0], [0.5, 0],
              [0, 0], [0, 0.5],
              [0, 0], [-0.5, 0]]
    assert dNtest == dNgood

    P = (1, 0)
    dNtest = geometry.SFQ8.dN(P)
    dNgood = [[0.5, 0], [-1, 0],
              [0.5, -0.5], [0.5, 0],
              [0.5, 0.5], [-1, 0],
              [0.5, 0], [-0.5, 0]]
    assert dNtest == dNgood

    P = (-1, 1)
    dNtest = geometry.SFQ8.dN(P)
    dNgood = [[0, 0.5], [0, 0],
              [0, 0], [0, 0],
              [-0.5, 0], [2, 0],
              [-1.5, 1.5], [0, -2]]
    assert dNtest == dNgood

    P = (0, 1)
    dNtest = geometry.SFQ8.dN(P)
    dNgood = [[0, 0.5], [0, -0.5],
              [0, 0.5], [0, -1],
              [0.5, 0.5], [0, 0.5],
              [-0.5, 0.5], [0, -1]]
    assert dNtest == dNgood

    P = (1, 1)
    dNtest = geometry.SFQ8.dN(P)
    dNgood = [[0, 0], [0, 0],
              [0, 0.5], [0, -2],
              [1.5, 1.5], [-2, 0],
              [0.5, 0], [0, 0]]
    assert dNtest == dNgood
