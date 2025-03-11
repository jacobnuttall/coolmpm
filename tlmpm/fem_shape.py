# MIT License
# 
# Copyright (c) 2025 Jacob Nuttall
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

# # # # # # # # # # # # # # # # # # # # # # # # # # 
# 
# 
#   CoolMPM : A Total Lagrangian MPM Code
# 
# 
# # # # # # # # # # # # # # # # # # # # # # # # # # 

import numpy as np
from functools import cache

@cache
def calc_GuassPoints(numips:int):
    """ Calculate the Guass-Legendre Quadrature points and weights 
    for a specified number of integration points.
    See: https://en.wikipedia.org/wiki/Gaussian_quadrature

    Args:
        numips (int): Number of integration points

    Returns:
        (points (np.array), weights (np.array)): the points and weights for numips
    """
    
    x = np.zeros(numips) # locations
    w = np.zeros(numips) # weights
    
    # Get Nth order Legendre Polynomial
    n = np.zeros(numips+1)
    n[-1] = 1
    Pn = np.polynomial.legendre.Legendre(n)
    x = Pn.roots()
    w = 2/(1-x**2)/(Pn.deriv()(x))**2
    
    return x, w

def XI(x,x0,x1):
    return 2*(x - 0.5*(x1 + x0))/(x1 - x0) 

def X(xi,x0,x1):
    return 0.5*(x1 - x0)*xi + 0.5*(x1 + x0)

def dXI_dx(x0,x1):
    return 2/(x1 - x0)

# def N(xi):
#    return 0.5*np.array([1-xi,1+xi])

# def B(xi):
#     return 0.5*np.array([-1.0,1.0])

def calc_element_nodes(n_en):
    return np.linspace(-1,1,n_en)

def ensure_dim_2d(x:np.array):
    x = np.asarray(x).astype(float)
    x = np.atleast_2d(x)
    if x.shape[1] != 1:
        x = x.T
    return x

def N(n_en, xi): 
    """ Calculate the finite element shape function

    Args:
        n_en (int): number of element nodes
        xi (np.ndarray | float): value(s) to evaluate at

    Returns:
        np.array: the finite element shape function
    """
    
    XI = calc_element_nodes(n_en)
    xi = ensure_dim_2d(xi)
    N = np.empty((n_en, xi.shape[0]))
    for i in range(n_en):
        denom = numer = 1
        for j in range(n_en):
            if j == i: continue
            denom *= XI[i] - XI[j]
            numer *= xi - XI[j]
        N[i] = np.squeeze(numer / denom)
    return np.array(N).T.squeeze()

def B(n_en, xi):
    """ Calculate the derivative of the finite element shape function

    Args:
        n_en (int): number of element nodes
        xi (np.ndarray | float): _description_

    Returns:
        np.array: the derivative of the finite element shape function
    """
    
    XI = calc_element_nodes(n_en)
    xi = ensure_dim_2d(xi)
    B = np.empty((n_en, xi.shape[0]))
    for i in range(n_en):
        numer = 0
        denom = 1
        for j in range(n_en):
            if j == i: continue
            numer_j = 1 + np.zeros_like(xi)
            denom *= XI[i] - XI[j]
            for k in range(n_en):
                if k == i or k == j: continue
                numer_j *= xi - XI[k]
                
            numer += numer_j
        B[i] = np.squeeze(numer/denom)

    return np.array(B).T.squeeze()