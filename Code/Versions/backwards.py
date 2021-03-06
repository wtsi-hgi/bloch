# Copyright (c) 2012 Genome Research Ltd. Author: Michelle Parker <mp18@sanger.ac.uk>
# This file is part of BLOCH.
# BLOCH is free software: you can redistribute it and/or modify it under the terms
# of the GNU General Public License as published by the Free Software Foundation;
# either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Generalquit() Public License for more
# details.
# You should have received a copy of the GNU General Public License along with
# this program.  If not, see <http://www.gnu.org/licenses/>.

#Import required python packages
from __future__ import division
import networkx as nx
import sys
import operator
import os
import math
import numpy as np
import itertools

from forward import *

#Function to find edge corresponding to matrix coordinate. g=flattened index. l=level(index of m)
def findedge(g,l):

    #Convert flattened index to coordinate
    t = np.unravel_index(g, (n[l],n[l]))

    f = [-1,-1]
    if l == 0:
        x = 1
    else:
        x = [j for j in range(glevel[l-1]+1,glevel[l]+1)]
        
    #Set f to equal tuple of ordered edges that the index g corresponds to
    for a, b in enumerate(G.out_edges(nbunch=x, keys=True, data=True)):
        if t[0] == a:
                f[0] = b
        if t[1] == a:
                f[1] = b

    if -1 in f:
        print 'Error in findedge'
    else:
        return f
   
for i in range(ll):
    print m[i]

#Initial probability
inprob = []

#Create probability distribution for sampling first allele
for i in m[3].flat:
    inprob.append(i/m[ll-1].sum())

#List of  flattened indices chosen from sampling 
s = [0]*ll

#Corresponding probabilities of the chosen positions 
p = [0]*ll

#Randomly choose first element in m[3] according to initial probabilities
s[ll-1] = np.random.choice(m[ll-1].size, p=inprob)

#Set probability of the chosen position
p[ll-1] = inprob[s[ll-1]]

#Iterate through levels in reverse order
for i in range(ll-1, 0, -1):
    prob = []

    #Set e to equal the edge description of edges sampled
    e = findedge(s[i], i)

    #For each position in forward probability matrix of level below, calculate sampling probabilities
    for b, j in enumerate(m[i-1].flat):
        d = findedge(b, i-1)
        prob.append((emission(GT[i],(e[0][3]['allele'],e[1][3]['allele']))*diptrans((e[0],d[0]),(e[1],d[1]))*j)/m[i].flat[s[i]])

    #Choose next sampled edge based on calculated probabilities
    s[i-1] = np.random.choice(m[i-1].size, p=prob)

    #Record probability of chosen edges.
    p[i-1] = prob[s[i-1]]

#Print descriptive list of chosen edges
for i in range(ll):
    print findedge(s[i], i)

print 'The pair of paths has sampling probability in either order of'

#Calculate probability of sampled path
print np.product(p)*2

