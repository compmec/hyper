# -*- coding: utf-8 -*-
#
# Script to test tensor module
#
"""
MNMNL Homework
Module 1: operations on tensors

Complete M1__tensor_template.py module and test your methods with this script.

You MUST verify your results either with a calculation "by hand" or with another
programm.
"""

#
# --- Namespace
#
# Import module tensor.py to access its methods with "tensor.method()".
# REQUIREMENTS:
# - the file M1__tensor_template.py must be renamed tensor.py,
# - it must be in the same folder as this script,
# - (or in a folder listed in the environment variable PYTHONPATH)
# WARNING:
# if you get an error such as "ImportError: No module named tensor", you have
# not fullfilled the requirements above.

import tensor

#
# --- Create a vector
#
v = tensor.vector()  # the only argument of the tensor.vector() method is d.
# In the definition of the function, the default value for d is 2 ("d=2").
# When the method is called in a script, one do not need to specify the value of d if one needs d=2.
v[0] = 1.0
v[1] = 2.0
# to access the elements of an array you must use brackets.
# Paranthesis are restricted to the call of functions.
print("v =\n", v)

#
# --- Test methods for 2nd-order tensors
#
print("\nTest methods for 2nd-order tensors")
print("==================================")
# Create a 2nd-order tensor, initialized to zero.
A = tensor.tensor()
print("\nA =\n", A)

# Compute identity 2nd-order tensor
I = tensor.I()
print("\nI =\n", I)

# Tensor product of two vectors
B = tensor.outerProd(v, v)
print("\nB =\n", B)

# Compute determinant
detB = tensor.det(B)
print("\ndetB =\n", detB)

# Compute trace
trB = tensor.trace(B)
print("\ntrB =\n", trB)

# Non symmetrical tensor
F = tensor.tensor()
F[0, 0] = 2.
F[0, 1] = 3.
F[1, 0] = 2
F[1, 1] = 5.
print("\nF = \n", F)
# Compute inverse
Finv = tensor.inv(F)
print("\nF^-1 =\n", Finv)

# Cauchy-Green tensors
C = tensor.rightCauchyGreen(F)  # C = F^t.F
print("\nC =\n", C)
b = tensor.leftCauchyGreen(F)   # b = F.F^t
print("\nb =\n", b)

#
# --- Test methods for 4th-order tensors
#
print("\nTest methods for 4th-order tensors")
print("==================================")
# Create a 4th-order tensor, initialized to zero
HH = tensor.tensor4()
print("\nHH =\n", HH)

# Compute identity 2nd-order tensor
II = tensor.II()
print("\nII =\n", II)

# Compute symmetrical identity 2nd-order tensor
# 0.5 * (delta_ik delta_jl + delta_il delta_jk)
IIs = tensor.IISym()
print("\nIIs =\n", IIs)

# Compute spherical operator
# delta_ij delta_kl
KK = tensor.KK()
print("\nKK =\n", KK)

# Compute 4th-order tensor product of two 2nd-order tensors
MM = tensor.outerProd4(C, C)
print("\nMM =\n", MM)
