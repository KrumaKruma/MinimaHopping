import numpy as np

class HistoryList:
    def __init__(self, ndim, nhist_max):
        self.ndim = ndim
        self.nhist_max = nhist_max
        self.histList = np.zeros((ndim, nhist_max))
        self.diffList = np.zeros((ndim, nhist_max))
        self.normalizedDiffList = np.zeros((ndim, nhist_max))
        self.icount = 0
        self.firstFull = True


    def add(self, vectorToAdd):
        if self.icount < self.nhist_max:
            self.histList[:, self.icount] = vectorToAdd
            self.icount += 1
            if self.icount > 1:
                self.diffList[:, self.icount - 2] = self.histList[:, self.icount - 1] - self.histList[:, self.icount - 2]
                self.normalizedDiffList[:, self.icount - 2] = self.diffList[:, self.icount - 2 ] / np.linalg.norm(self.diffList[:, self.icount - 2 ])
            return self.icount - 1
        else:
            # print('list full')
            # make place in history list
            self.histList[:, :(self.nhist_max-1)] = self.histList[:, 1:(self.nhist_max)]
            # add last element
            self.histList[:, -1] = vectorToAdd

            if self.firstFull:
                self.firstFull = False
                self.diffList[:, -1] = self.diffList[:, -1] = self.histList[:, -1] - self.histList[:,-2]
                self.normalizedDiffList[:, -1] = self.diffList[:, -1] / np.linalg.norm(self.diffList[:, -1])
            else: 
                self.diffList[:, :(self.nhist_max-1)] = self.diffList[:, 1:(self.nhist_max)]
                self.diffList[:, -1] = self.histList[:, -1] - self.histList[:,-2]
                self.normalizedDiffList[:, :(self.nhist_max-1)] = self.normalizedDiffList[:, 1:(self.nhist_max)]
                self.normalizedDiffList[:, -1] = self.diffList[:, -1] / np.linalg.norm(self.diffList[:, -1])
            return self.nhist_max

    def print_me(self, n):
        print('histlist')
        print(self.histList)
        print('difflist')
        print(self.diffList[:, :(n)])
        print('normalizeddifflist')
        print(self.normalizedDiffList[:, :(n)])


#a = HistoryList(2, 5)
#for i in range(1, 10):
#    x = np.ones(2) * i**2
#    n = a.add(x)
#    a.print_me(n)
#    print('length', n)
#    print('')