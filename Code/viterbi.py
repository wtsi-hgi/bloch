# Copyright (c) 2012 Genome Research Ltd. Author: Michelle Parker <mp18@sanger.ac.uk>
# This file is part of BLOCH.
# BLOCH is free software: you can redistribute it and/or modify it under the terms
# of the GNU General Public License as published by the Free Software Foundation;
# either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
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

#v is the list of matrices of viterbi probabilities at each level
v = [np.zeros(shape=(n[0],n[0]))]

#arglist is the list of matrices of path to find maximum viterbi probabilities
arglist = [np.zeros((n[0],n[0],2))]
arglist[0].fill(-1)

#Append matrices to list arglist
for i in range(ll-1):
    v.append(np.zeros(shape=(n[i+1],n[i+1])))
    arglist.append(np.zeros((n[i+1],n[i+1],2)))
    arglist[i+1].fill(-1)


#Initiation. Iterate through pairs of outgoing edges from node 1.
for a, b  in itertools.product([(i, j) for i, j in enumerate(G.out_edges(1, keys=True, data=True))], repeat=2):
    
    #If emmsion probability does not equal 0
    if emission(GT[0],(a[1][3]['allele'],b[1][3]['allele'])) != 0:        
        
        if a[1] == b[1]:
            t = hapinitial(a[1][3]['allele'])
            var = t*t            
        else:
            var = dipinitial(a[1][3]['allele'], b[1][3]['allele'])

        #Matrix element is set to calculated diploid initial probability
        v[0][a[0]][b[0]] = var

#Induction. Iterate through each level
for i in range(1,ll):
    
    #Iterate through ordered pairs of outgoing edges in each level
    for a, b in itertools.product([(c, d) for c, d in enumerate(G.out_edges(nbunch=[j for j in range(glevel[i-1]+1,glevel[i]+1)], keys=True, data=True))], repeat=2):

        #If emission probability does not equal 0
        if emission(GT[i],(a[1][3]['allele'],b[1][3]['allele'])) != 0:
            
            if i == 1:
                nodes = 1
            else:
                nodes = [j for j in range(glevel[i-2]+1,glevel[i-1]+1)]

            max = 0
            args = [-1., -1.]
            for c, d  in itertools.product([(e, f) for e, f in enumerate(G.out_edges(nbunch=nodes, keys=True, data=True))], repeat=2):

                value = (v[i-1][c[0]][d[0]]*diptrans([a[1],c[1]], [b[1],d[1]]))
               
                if max < value:
                    max = value
                    args = [c[0],d[0]]
                    
            #Matrix element is set to var calculation formula
            v[i][a[0]][b[0]] = max
            arglist[i][a[0]][b[0]] = args

for i in range(ll):
    print 'level '+str(i)
    print m[i]
    print v[i]
    print arglist[i]
