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

class HistoryList:
    """Historylist that is used by the sqnm class.
    More informations about the SQNM algorithm can be found here: https://aip.scitation.org/doi/10.1063/1.4905665
    """

    def __init__(self, ndim, nhist_max):
        self.ndim = ndim
        self.nhist_max = nhist_max
        self.histList = np.zeros((ndim, nhist_max))
        self.diffList = np.zeros((ndim, nhist_max))
        self.normalizedDiffList = np.zeros((ndim, nhist_max))
        self.icount = 0
        self.firstFull = True


    def add(self, vectorToAdd):
        """
        Add a vector to the history list.
        """

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

    def _print_me(self, n):
        print('histlist')
        print(self.histList)
        print('difflist')
        print(self.diffList[:, :(n)])
        print('normalizeddifflist')
        print(self.normalizedDiffList[:, :(n)])
