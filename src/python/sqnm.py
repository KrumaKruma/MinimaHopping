#!/usr/bin/env python3

# The variable cell shape optimization method is based on the following 
# paper: https://arxiv.org/abs/2206.07339
# More details about the SQNM optimization method are available here:
# https://comphys.unibas.ch/publications/Schaefer2015.pdf
# Author of this document: Moritz Gubler 

# Copyright (C) 2022 Moritz Gubler
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import historylist
import time

class SQNM:
    """This class is an implementation of the stabilized quasi newton optimization method.
    More informations about the algorithm can be found here: https://aip.scitation.org/doi/10.1063/1.4905665
    """
    
    def __init__(self, ndim, nhist_max, alpha, eps_supsp, alpha_min):
        """
        Parameters
        ----------
        ndim : int
            Number of dimensions of the optimization problem
        nhist_max : int
            maximal length of history list
        alpha: double
            Initial step size. Should be approximately the inverse of the largest eigenvalue of the Hessian matrix.
            If alpha is negative, the inial step size is estimated using the mechanism from section 6.4 of the 
            vc-sqnm paper: https://arxiv.org/abs/2206.07339
            beta will then be equal to minus alpha. Good choices for beta are 0.1 in hartee / bohr^2 and
            0.001 in eV / A^2
        eps_subsp: double
            Lower limit on linear dependencies in history list.
        alpha_min: double
            Lowest step size that is allowed.
        """

        self.ndim = ndim
        self.nhist_max = nhist_max
        if ndim < nhist_max:
            print('ndim > nhist_max. Seting nhist_max to ndim')
            self.nhist_max = self.ndim
        self.eps_subsp = eps_supsp
        self.alpha_min = alpha_min
        self.xlist = historylist.HistoryList(self.ndim, self.nhist_max)
        self.flist = historylist.HistoryList(self.ndim, self.nhist_max)
        
        # when alpha is smaller than zero estimate initial step size.
        self.estimate_step_size = False
        if alpha <= 0:
            self.estimate_step_size = True
            self.alpha = - alpha
        else:
            self.alpha = alpha
        
        self.dir_of_descent = np.zeros(ndim)
        self.expected_positions = np.zeros(ndim)
        self.prev_f_of_x = 0.0
        self.prev_df_dx = np.zeros(ndim)
        self.s_evec = np.zeros((nhist_max, nhist_max))
        self.s_eval = np.zeros(nhist_max)
        self.dr_subsp = np.zeros((ndim, nhist_max))
        self.df_subsp = np.zeros((ndim, nhist_max))
        self.h_evec_subsp = np.zeros((nhist_max, nhist_max))
        self.h_evec = np.zeros((ndim, nhist_max))
        self.h_eval = np.zeros(nhist_max)
        self.res = np.zeros(nhist_max)
        self.gainratio = 0.0
        self.nhist = 0

    def sqnm_step(self, x, f_of_x, df_dx):
        """ Calculates a set of new coordinates based on the function value and derivatives provide on input.

        The idea behind this function is, that the user evaluates the function at the new point this method suggested and
        then calls this method again with the function value at the new point until convergence was reached.

        Parameters
        ----------
        x: numpy array
            Numpy array containg position vector x.
        f_of_x: double
            Value of target function at point x.
        df_dx: numpy array
            derivative of function f with respect to x.
        """

        # check if norm of gradient is already zero:
        # and return zero if this is the case
        if np.linalg.norm(df_dx) < 10e-13:
            self.dir_of_descent = np.zeros(self.ndim)
            return self.dir_of_descent

        self.nhist = self.xlist.add(x)
        self.flist.add(df_dx)

        # check if first step
        if self.nhist == 0:
            self.dir_of_descent = -self.alpha * df_dx
        else:

            # check if positions have been modified
            if np.max(np.abs(x - self.expected_positions)) > 10.0e-9:
                print("SQNM was not called with positions that were expected. If this was not done on purpose, it is probably a bug.")
                print("Were atoms that left the simulation box put back into the cell? This is not allowed.")

            if self.estimate_step_size:
                l1 = (f_of_x - self.prev_f_of_x + self.alpha * np.linalg.norm(self.prev_df_dx)**2) / (0.5 * (self.alpha**2) * np.linalg.norm(self.prev_df_dx)**2)
                l2 = np.linalg.norm(df_dx - self.prev_df_dx) / (self.alpha * np.linalg.norm(self.prev_df_dx))
                self.alpha = min(1/ l1, 1/l2)
                print("      Automatic initial step size guess for geometry optimization: ", self.alpha)
                self.estimate_step_size = False
            else:
                # calculate and adjust gainratio
                self.gainratio = (f_of_x - self.prev_f_of_x) / ( .5 * np.dot(self.dir_of_descent, self.prev_df_dx) )
                if not self.estimate_step_size:
                    if self.gainratio < .5:
                        self.alpha = max(self.alpha_min, self.alpha * .65)
                    if self.gainratio > 1.05:
                        self.alpha *= 1.05

            # calculate overlap matrix of basis
            self.s_evec[:self.nhist, :self.nhist] = self.xlist.normalizedDiffList[:, :self.nhist].T \
                @ self.xlist.normalizedDiffList[:, :self.nhist]
            self.s_eval[:self.nhist], self.s_evec[:self.nhist, :self.nhist] \
                = np.linalg.eigh(self.s_evec[:self.nhist, :self.nhist])

            # remove noisy directions from subspace
            dim_subsp = sum(self.s_eval[:self.nhist] / self.s_eval[self.nhist - 1] > self.eps_subsp)
            #dim_subsp = min(1, dim_subsp)
            self.s_eval[:dim_subsp] = self.s_eval[(self.nhist - dim_subsp):self.nhist]
            self.s_evec[:, :dim_subsp] = self.s_evec[:, (self.nhist - dim_subsp):self.nhist]

            # compute eq. 11
            self.dr_subsp[:, :dim_subsp] = np.einsum('hi,kh->ki', self.s_evec[:self.nhist,:dim_subsp], self.xlist.normalizedDiffList[:, :self.nhist])\
                / np.sqrt(self.s_eval[:dim_subsp])
            self.df_subsp[:, :dim_subsp] = np.einsum('hi,kh,h->ki', self.s_evec[:self.nhist,:dim_subsp], self.flist.diffList[:, :self.nhist], np.divide(1.0, np.linalg.norm(self.xlist.diffList[:, :self.nhist], axis=0)) )\
                / np.sqrt(self.s_eval[:dim_subsp])
                    
            # compute eq 13
            self.h_evec_subsp[:dim_subsp, :dim_subsp] = .5 * (self.df_subsp[:, :dim_subsp].T @ self.dr_subsp[:, :dim_subsp] \
                + self.dr_subsp[:, :dim_subsp].T @ self.df_subsp[:, :dim_subsp] )
            self.h_eval[:dim_subsp], self.h_evec_subsp[:dim_subsp, :dim_subsp] \
                = np.linalg.eigh(self.h_evec_subsp[:dim_subsp, :dim_subsp])
            
            # compute eq. 15
            self.h_evec[:, :dim_subsp] = np.einsum('ki,hk->hi', self.h_evec_subsp[:dim_subsp, :dim_subsp], self.dr_subsp[:, :dim_subsp])

            # compute eq. 20
            self.res[:dim_subsp] = np.linalg.norm(
                - self.h_eval[:dim_subsp] * self.h_evec[:, :dim_subsp] 
                + self.df_subsp[:, :dim_subsp] @ self.h_evec_subsp[:dim_subsp, :dim_subsp]
                , axis = 0)

            # modify eigenvalues according to eq. 18
            self.h_eval[:dim_subsp] = np.sqrt(self.h_eval[:dim_subsp]**2 + self.res[:dim_subsp]**2)

            # decompose gradient according to eq. 16
            self.dir_of_descent = (df_dx - np.einsum('i, ki -> k', self.h_evec[:, :dim_subsp].T @ df_dx, self.h_evec[:, :dim_subsp])) * self.alpha

            # apply preconditioning to subspace gradient (eq. 21)
            self.dir_of_descent = self.dir_of_descent + \
                np.einsum('i, ki, i -> k', self.h_evec[:,:dim_subsp].T @ df_dx\
                , self.h_evec[:, :dim_subsp], np.divide(1.0, self.h_eval[:dim_subsp]) )

            self.dir_of_descent = - self.dir_of_descent
        self.expected_positions = x + self.dir_of_descent
        self.prev_f_of_x = f_of_x
        self.prev_df_dx = df_dx
        return self.dir_of_descent

    def lower_bound(self):
        """ Returns an estimate of a lower bound for the local minumum.
        The estimate is only accurate when the optimization is converged.
        """

        if self.nhist < 1:
            print("At least one optimization step needs to be done before lower_limit can be called.")
            return 0
        return self.prev_f_of_x - .5 * np.dot(self.prev_df_dx, self.prev_df_dx) / self.h_eval[0]

# the rest of this file can be used for testing only

def _test_fun(x):
    n = len(x)
    a = np.zeros((n, n))
    for i in range(n):
        a[i, i] = .5*i + 1
    f = .5 * x.T @ a @ x
    df = a @ x
    return f, df


def _tests():
    n = 400
    nhistx = 10
    alpha = .005
    x = np.zeros(n)
    x[:] = 1.0
    x[1] = 1.4

    opt = SQNM(n, nhistx, alpha, 1e-4, 1e-2)


    for i in range(100):
        f, df = _test_fun(x)
        #print(i)
        #print('norm of x', np.linalg.norm(x))
        #print('f(x)', f)
        #print('norm of forces', np.linalg.norm(df))
        #print('')
        t1 = time.time()
        dd = opt.sqnm_step(x, f, df)
        t2 = time.time()
        print(f, opt.lower_limit()) 
        x += dd


if __name__ == "__main__":
    _tests()