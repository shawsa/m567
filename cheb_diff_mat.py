'''
Original comments in MATLAB
%CHEBDIFFMAT  Chebyshev pseudospectral differentiation matrices over [a,b]
%  The function [x, DM] =  chebdiffmat(N,M,a,b) computes the differentiation
%  matrices D1, D2, ..., DM on Chebyshev nodes spaced over the interval
%  [a,b].
%
%  Input:
%  N:        Size of differentiation matrix.
%  M:        Number of derivatives required (integer).
%  [a,b]:    Interval for the approximation.
%  Note:     0 < M <= N-1.
%
%  Output:
%  DM:       DM(1:N,1:N,ell) contains ell-th derivative matrix, ell=1..M.
%
%  The code implements two strategies for enhanced accuracy suggested by
%          W.Don and S. Solomonoff in SIAM J. Sci. Comp. Vol. 6, pp.
%          1253--1268 (1994).
%  The two strategies are (a) the use of trigonometric identities
%  to avoid the computation of differences x(k)-x(j) and (b) the use of the
%  "flipping trick" which is necessary since sin t can be computed to high
%  relative precision when t is small whereas sin (pi-t) cannot.
%
%  Note added May 2003:  It may, in fact, be slightly better not to
%  implement the strategies (a) and (b).   Please consult the following
%  paper for details:
%          R. Baltensperger and M.R. Trummer. Spectral Differencing with a
%          Twist, by , to appear in SIAM J. Sci. Comp.
%
%  Original author: J.A.C. Weideman, S.C. Reddy 1998.  Help notes modified
%  by JACW, May 2003.
%
%  Modified by G.B. Wright 2018 to scale to any interval, and a few minor
%  improvements for later versions of MATLAB.

Converted to Python3+Numpy by S. Shaw 2018



'''


import numpy as np

def cheb_diff_mat(n, m, a, b):
    thetas = np.linspace(0, np.pi, n, endpoint=True)
    T = np.repeat(.5*thetas.reshape((n,1)), n, axis=1)
    DX = 2 * np.sin(T.T + T) * np.sin(T.T - T)
    DX[n//2:] = -np.rot90(DX[:(n+1)//2], 2)
    DX[np.diag_indices(n)] = np.ones(n)
    C = (-1.0) ** np.add.outer(np.arange(n), np.arange(n))
    C[0]    *= 2
    C[-1]   *= 2
    C[:,0]  *= .5
    C[:,-1] *= .5
    Z = 1/DX
    Z[np.diag_indices(n)] = np.zeros(n)
    D = np.eye(n)
    DM = np.zeros((m,n,n))
    scale_factor = 1
    for i in range(m):
        D = (i+1)*Z*(C * np.multiply.outer(np.diag(D), np.ones(n)) - D)
        D[np.diag_indices(n)] = -D @ np.ones(n)
        scale_factor *= 2/(b-a)
        DM[i] = scale_factor * D[::-1,::-1] #account for flipping the points
    xs = -np.cos(thetas) #standard Chebyshev points
    xs = (xs+1) * (b-a)/2 + a #shifted
    return xs, DM
