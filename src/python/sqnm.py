import numpy as np
import historylist

class SQNM:
    def __init__(self, ndim, nhist_max, alpha, eps_supsp, alpha_min):
        self.ndim = ndim
        self.nhist_max = nhist_max
        self.eps_subsp = eps_supsp
        self.alpha_min = alpha_min
        self.xlist = historylist.HistoryList(self.ndim, self.nhist_max)
        self.flist = historylist.HistoryList(self.ndim, self.nhist_max)
        self.alpha = alpha
        self.dir_of_descent = np.zeros(ndim)
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
        self.res_temp = np.zeros(ndim)
        self.gainratio = 0.0
        self.nhist = 0

    def sqnm_step(self, x, f_of_x, df_dx):
        self.nhist = self.xlist.add(x)
        self.flist.add(df_dx)

        # check if first step
        if self.nhist == 0:
            self.dir_of_descent = -self.alpha * df_dx
        else:
            # calculate and adjust gainratio
            self.gainratio = (f_of_x - self.prev_f_of_x) / ( .5 * np.dot(self.dir_of_descent, self.prev_df_dx) )
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
            dim_subsp = 0
            for i in range(self.nhist):
                if self.s_eval[i] / self.s_eval[self.nhist - 1] > self.eps_subsp:
                    dim_subsp += 1

            self.s_eval[:dim_subsp] = self.s_eval[(self.nhist - dim_subsp):self.nhist]
            self.s_evec[:, :dim_subsp] = self.s_evec[:, (self.nhist - dim_subsp):self.nhist]

            # compute eq. 11
            self.dr_subsp[:, :] = 0.0
            self.df_subsp[:, :] = 0.0
            for i in range(dim_subsp):
                for ihist in range(self.nhist):
                    self.dr_subsp[:, i] += self.s_evec[ihist, i] * self.xlist.normalizedDiffList[:, ihist]
                    self.df_subsp[:, i] += self.s_evec[ihist, i] * self.flist.diffList[:, ihist] \
                        / np.linalg.norm(self.xlist.diffList[:, ihist])
                self.dr_subsp[:, i] = self.dr_subsp[:, i] / np.sqrt(self.s_eval[i])
                self.df_subsp[:, i] = self.df_subsp[:, i] / np.sqrt(self.s_eval[i])
                    
            # compute eq 13
            self.h_evec_subsp[:dim_subsp, :dim_subsp] = .5 * (self.df_subsp[:, :dim_subsp].T @ self.dr_subsp[:, :dim_subsp] \
                + self.dr_subsp[:, :dim_subsp].T @ self.df_subsp[:, :dim_subsp] )
            self.h_eval[:dim_subsp], self.h_evec_subsp[:dim_subsp, :dim_subsp] \
                = np.linalg.eigh(self.h_evec_subsp[:dim_subsp, :dim_subsp])
            
            # compute eq. 15
            self.h_evec[:, :] = 0.0
            for i in range(dim_subsp):
                for k in range(dim_subsp):
                    self.h_evec[:, i] += self.h_evec_subsp[k, i] * self.dr_subsp[:, k]

            # compute eq. 20
            for j in range(dim_subsp):
                self.res_temp = - self.h_eval[j] * self.h_evec[:, j]
                for k in range(dim_subsp):
                    self.res_temp += self.h_evec_subsp[k, j] * self.df_subsp[:, k]
                self.res[j] = np.linalg.norm(self.res_temp)

            # modify eigenvalues according to eq. 18
            self.h_eval[:dim_subsp] = np.sqrt(self.h_eval[:dim_subsp]**2 + self.res[:dim_subsp]**2)

            # decompose gradient according to eq. 16
            self.dir_of_descent = df_dx
            for i in range(dim_subsp):
                self.dir_of_descent = self.dir_of_descent - np.dot(self.h_evec[:, i], df_dx) * self.h_evec[:, i]
            self.dir_of_descent = self.dir_of_descent * self.alpha

            # apply preconditioning to subspace gradient (eq. 21)
            for idim in range(dim_subsp):
                self.dir_of_descent = self.dir_of_descent + np.dot( df_dx, self.h_evec[:, idim] ) * self.h_evec[:, idim] / self.h_eval[idim]
            #print('dd2', self.dir_of_descent)

            self.dir_of_descent = - self.dir_of_descent
        self.prev_f_of_x = f_of_x
        self.prev_df_dx = df_dx
        return self.dir_of_descent


"""
def test_fun(x):
    n = len(x)
    a = np.zeros((n, n))
    for i in range(n):
        a[i, i] = .5*i + 1
    f = .5 * x.T @ a @ x
    df = a @ x
    return f, df



n = 400
nhistx = 20
alpha = .1
x = np.zeros(n)
x[:] = 1.0
x[1] = 1.4

opt = SQNM(n, nhistx, alpha, 1e-4, 1e-2)


for i in range(40):
    f, df = test_fun(x)
    print(i)
    print('norm of x', np.linalg.norm(x))
    print(f)
    print('norm of forces', np.linalg.norm(df))
    print('')
    t1 = time.time()
    dd = opt.sqnm_step(x, f, df)
    t2 = time.time()
    print('time', t2 -t1)
    # print('dd', dd)
    x += dd
"""