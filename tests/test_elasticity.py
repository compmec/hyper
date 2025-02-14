# , -*- coding: utf-8, -*-
"""
"""

import numpy as np
import pytest

from hyper import elasticity, tensor


def test_StVenantKirchhoffValues():
    E = 210
    nu = 0.3
    G = 1050 / 13
    mu = 1050 / 13
    lamda = 1575 / 13
    material = elasticity.StVenantKirchhoffElasticity(E, nu)
    assert material.E == pytest.approx(E)
    assert material.nu == pytest.approx(nu)
    assert material.G == pytest.approx(G)
    assert material.LAME1 == pytest.approx(lamda)
    assert material.LAME2 == pytest.approx(mu)

    E = 26
    nu = 0.3
    G = 10
    mu = 10
    lamda = 15
    material = elasticity.StVenantKirchhoffElasticity(E, nu)
    assert material.E == pytest.approx(E)
    assert material.nu == pytest.approx(nu)
    assert material.G == pytest.approx(G)
    assert material.LAME1 == pytest.approx(lamda)
    assert material.LAME2 == pytest.approx(mu)

    E = 42
    nu = 0.4
    G = 15
    mu = 15
    lamda = 60
    material = elasticity.StVenantKirchhoffElasticity(E, nu)
    assert material.E == pytest.approx(E)
    assert material.nu == pytest.approx(nu)
    assert material.G == pytest.approx(G)
    assert material.LAME1 == pytest.approx(lamda)
    assert material.LAME2 == pytest.approx(mu)


def test_StVenantKirchhoffPotential():
    E = 52
    nu = 0.3
    material = elasticity.StVenantKirchhoffElasticity(E, nu)

    C = np.eye(2)  # No deformation at all
    phitest = material.potential(C)
    assert phitest == 0

    C = np.eye(3)  # No deformation at all
    phitest = material.potential(C)
    assert phitest == 0

    C = np.ones((2, 2))
    phitest = material.potential(C)
    assert phitest == 10

    C = 2 * np.eye(2) + 1
    phitest = material.potential(C)
    assert phitest == 110

    C = np.eye(2) + 1
    phitest = material.potential(C)
    assert phitest == 35


def test_StVenantKirchhoffStress():
    E = 26
    nu = 0.3
    material = elasticity.StVenantKirchhoffElasticity(E, nu)

    C = np.eye(2)  # No deformation at all
    PK2test = material.stress(C)
    assert PK2test.ndim == 2
    assert PK2test.shape == (2, 2)
    np.testing.assert_almost_equal(PK2test, np.zeros((2, 2)))

    C = np.eye(3)  # No deformation at all
    PK2test = material.stress(C)
    assert PK2test.ndim == 2
    assert PK2test.shape == (3, 3)
    np.testing.assert_almost_equal(PK2test, np.zeros((3, 3)))

    C = np.ones((2, 2))
    PK2test = material.stress(C)
    PK2good = 10 * np.array([[0, 1], [1, 0]])
    assert PK2test.ndim == 2
    assert PK2test.shape == (2, 2)
    np.testing.assert_almost_equal(PK2test, PK2good)

    C = 2 * np.eye(2) + 1
    PK2test = material.stress(C)
    PK2good = 10 * np.array([[5, 1], [1, 5]])
    assert PK2test.ndim == 2
    assert PK2test.shape == (2, 2)
    np.testing.assert_almost_equal(PK2test, PK2good)

    C = np.eye(2) + 1
    PK2test = material.stress(C)
    PK2good = 5 * np.array([[5, 2], [2, 5]])
    assert PK2test.ndim == 2
    assert PK2test.shape == (2, 2)
    np.testing.assert_almost_equal(PK2test, PK2good)


def test_StVenantKirchhoffStiffness():
    E = 26
    nu = 0.3
    material = elasticity.StVenantKirchhoffElasticity(E, nu)

    Mgood = [
        [[[7, 0], [0, 3]], [[0, 2], [2, 0]]],
        [[[0, 2], [2, 0]], [[3, 0], [0, 7]]],
    ]
    Mgood = 5 * np.array(Mgood)

    C = np.eye(2)  # No deformation at all
    Mtest = material.stiffness(C)
    assert Mtest.ndim == 4
    assert Mtest.shape == (2, 2, 2, 2)
    np.testing.assert_almost_equal(Mtest, Mgood)

    C = np.ones((2, 2))
    Mtest = material.stiffness(C)
    assert Mtest.ndim == 4
    assert Mtest.shape == (2, 2, 2, 2)
    np.testing.assert_almost_equal(Mtest, Mgood)

    C = 2 * np.eye(2) + 1
    Mtest = material.stiffness(C)
    assert Mtest.ndim == 4
    assert Mtest.shape == (2, 2, 2, 2)
    np.testing.assert_almost_equal(Mtest, Mgood)

    C = np.eye(2) + 1
    Mtest = material.stiffness(C)
    assert Mtest.ndim == 4
    assert Mtest.shape == (2, 2, 2, 2)
    np.testing.assert_almost_equal(Mtest, Mgood)

    C = np.eye(3)  # No deformation at all
    Mgood = [
        [
            [[7, 0, 0], [0, 3, 0], [0, 0, 3]],
            [[0, 2, 0], [2, 0, 0], [0, 0, 0]],
            [[0, 0, 2], [0, 0, 0], [2, 0, 0]],
        ],
        [
            [[0, 2, 0], [2, 0, 0], [0, 0, 0]],
            [[3, 0, 0], [0, 7, 0], [0, 0, 3]],
            [[0, 0, 0], [0, 0, 2], [0, 2, 0]],
        ],
        [
            [[0, 0, 2], [0, 0, 0], [2, 0, 0]],
            [[0, 0, 0], [0, 0, 2], [0, 2, 0]],
            [[3, 0, 0], [0, 3, 0], [0, 0, 7]],
        ],
    ]
    Mgood = 5 * np.array(Mgood)
    Mtest = material.stiffness(C)
    assert Mtest.ndim == 4
    assert Mtest.shape == (3, 3, 3, 3)
    np.testing.assert_almost_equal(Mtest, Mgood)


def test_NeoHookeanValues():
    E = 210
    nu = 0.3
    G = 1050 / 13
    mu = 1050 / 13
    lamda = 1575 / 13
    material = elasticity.NeoHookeanElasticity(E, nu)
    assert material.E == pytest.approx(E)
    assert material.nu == pytest.approx(nu)
    assert material.G == pytest.approx(G)
    assert material.LAME1 == pytest.approx(lamda)
    assert material.LAME2 == pytest.approx(mu)

    E = 26
    nu = 0.3
    G = 10
    mu = 10
    lamda = 15
    material = elasticity.NeoHookeanElasticity(E, nu)
    assert material.E == pytest.approx(E)
    assert material.nu == pytest.approx(nu)
    assert material.G == pytest.approx(G)
    assert material.LAME1 == pytest.approx(lamda)
    assert material.LAME2 == pytest.approx(mu)

    E = 42
    nu = 0.4
    G = 15
    mu = 15
    lamda = 60
    material = elasticity.NeoHookeanElasticity(E, nu)
    assert material.E == pytest.approx(E)
    assert material.nu == pytest.approx(nu)
    assert material.G == pytest.approx(G)
    assert material.LAME1 == pytest.approx(lamda)
    assert material.LAME2 == pytest.approx(mu)


def test_NeoHookeanPotential():
    E = 52
    nu = 0.3
    material = elasticity.NeoHookeanElasticity(E, nu)

    C = np.eye(3)  # No deformation at all
    phitest = material.potential(C)
    assert phitest == pytest.approx(0)

    C = np.eye(3)  # No deformation at all
    phitest = material.potential(C)
    assert phitest == pytest.approx(0)

    C = 2 * np.eye(3) + 1
    phitest = material.potential(C)
    assert phitest == pytest.approx(63.6967217)

    C = np.eye(3) + 1
    phitest = material.potential(C)
    assert phitest == pytest.approx(23.3438516)


def test_NeoHookeanStress():
    E = 26
    nu = 0.3
    material = elasticity.NeoHookeanElasticity(E, nu)

    C = np.eye(2)  # No deformation at all
    PK2test = material.stress(C)
    PK2good = np.zeros((2, 2))
    assert PK2test.ndim == 2
    assert PK2test.shape == (2, 2)
    np.testing.assert_almost_equal(PK2test, PK2good)

    C = np.eye(3)  # No deformation at all
    PK2test = material.stress(C)
    PK2good = np.zeros((3, 3))
    assert PK2test.ndim == 2
    assert PK2test.shape == (3, 3)
    np.testing.assert_almost_equal(PK2test, PK2good)

    C = 0.1 * np.eye(2) + 1
    PK2test = material.stress(C)
    PK2good = [[-103.692, 103.356], [103.356, -103.692]]
    assert PK2test.ndim == 2
    assert PK2test.shape == (2, 2)
    np.testing.assert_almost_equal(PK2test, PK2good, decimal=3)

    C = 2 * np.eye(2) + 1
    PK2test = material.stress(C)
    PK2good = [[12.098, -0.699], [-0.699, 12.098]]
    assert PK2test.ndim == 2
    assert PK2test.shape == (2, 2)
    np.testing.assert_almost_equal(PK2test, PK2good, decimal=3)

    C = np.eye(2) + 1
    PK2test = material.stress(C)
    PK2good = [[8.826, 0.586], [0.586, 8.826]]
    assert PK2test.ndim == 2
    assert PK2test.shape == (2, 2)
    np.testing.assert_almost_equal(PK2test, PK2good, decimal=3)


def test_NeoHookeanStiffness():
    E = 26
    nu = 0.3
    material = elasticity.NeoHookeanElasticity(E, nu)

    C = np.eye(2)  # No deformation at all
    Mgood = [
        [[[7, 0], [0, 3]], [[0, 2], [2, 0]]],
        [[[0, 2], [2, 0]], [[3, 0], [0, 7]]],
    ]
    Mgood = 5 * np.array(Mgood)
    Mtest = material.stiffness(C)
    assert Mtest.ndim == 4
    assert Mtest.shape == (2, 2, 2, 2)
    np.testing.assert_almost_equal(Mtest, Mgood)

    C = 0.1 * np.eye(2) + 1
    Mgood = [
        [
            [[1602.624, -1456.93], [-1456.93, 1395.911]],
            [[-1456.93, 1427.839], [1427.839, -1456.93]],
        ],
        [
            [[-1456.93, 1427.839], [1427.839, -1456.93]],
            [[1395.911, -1456.93], [-1456.93, 1602.624]],
        ],
    ]
    Mtest = material.stiffness(C)
    assert Mtest.ndim == 4
    assert Mtest.shape == (2, 2, 2, 2)
    np.testing.assert_almost_equal(Mtest, Mgood, decimal=2)

    C = 2 * np.eye(2) + 1
    Mgood = [
        [
            [[0.535, -0.178], [-0.178, 1.934]],
            [[-0.178, -0.639], [-0.639, -0.178]],
        ],
        [
            [[-0.178, -0.639], [-0.639, -0.178]],
            [[1.934, -0.178], [-0.178, 0.535]],
        ],
    ]
    Mtest = material.stiffness(C)
    assert Mtest.ndim == 4
    assert Mtest.shape == (2, 2, 2, 2)
    np.testing.assert_almost_equal(Mtest, Mgood, decimal=3)

    C = np.eye(2) + 1
    Mgood = [
        [
            [[8.231, -4.115], [-4.115, 7.057]],
            [[-4.115, 2.644], [2.644, -4.115]],
        ],
        [
            [[-4.115, 2.644], [2.644, -4.115]],
            [[7.057, -4.115], [-4.115, 8.231]],
        ],
    ]
    Mtest = material.stiffness(C)
    assert Mtest.ndim == 4
    assert Mtest.shape == (2, 2, 2, 2)
    np.testing.assert_almost_equal(Mtest, Mgood, decimal=3)

    C = np.eye(3)  # No deformation at all
    Mgood = [
        [
            [[7, 0, 0], [0, 3, 0], [0, 0, 3]],
            [[0, 2, 0], [2, 0, 0], [0, 0, 0]],
            [[0, 0, 2], [0, 0, 0], [2, 0, 0]],
        ],
        [
            [[0, 2, 0], [2, 0, 0], [0, 0, 0]],
            [[3, 0, 0], [0, 7, 0], [0, 0, 3]],
            [[0, 0, 0], [0, 0, 2], [0, 2, 0]],
        ],
        [
            [[0, 0, 2], [0, 0, 0], [2, 0, 0]],
            [[0, 0, 0], [0, 0, 2], [0, 2, 0]],
            [[3, 0, 0], [0, 3, 0], [0, 0, 7]],
        ],
    ]
    Mgood = 5 * np.array(Mgood)
    Mtest = material.stiffness(C)
    assert Mtest.ndim == 4
    assert Mtest.shape == (3, 3, 3, 3)
    np.testing.assert_almost_equal(Mtest, Mgood)

    C = 2 * np.eye(3) + 0.5
    Mgood = [
        [
            [[-0.84, 0.14, 0.14], [0.14, 2.65, -0.55], [0.14, -0.55, 2.65]],
            [[0.14, -1.77, 0.33], [-1.77, 0.14, 0.33], [0.33, 0.33, -0.56]],
            [[0.14, 0.33, -1.77], [0.33, -0.56, 0.33], [-1.77, 0.33, 0.14]],
        ],
        [
            [[0.14, -1.77, 0.33], [-1.77, 0.14, 0.33], [0.33, 0.33, -0.56]],
            [[2.66, 0.14, -0.56], [0.14, -0.84, 0.14], [-0.56, 0.14, 2.66]],
            [[-0.56, 0.33, 0.33], [0.33, 0.14, -1.77], [0.33, -1.77, 0.14]],
        ],
        [
            [[0.14, 0.33, -1.77], [0.33, -0.56, 0.33], [-1.77, 0.33, 0.14]],
            [[-0.56, 0.33, 0.33], [0.33, 0.14, -1.77], [0.33, -1.77, 0.14]],
            [[2.66, -0.56, 0.14], [-0.56, 2.66, 0.14], [0.14, 0.14, -0.84]],
        ],
    ]
    Mtest = material.stiffness(C)
    assert Mtest.ndim == 4
    assert Mtest.shape == (3, 3, 3, 3)
    np.testing.assert_almost_equal(Mtest, Mgood, decimal=2)
