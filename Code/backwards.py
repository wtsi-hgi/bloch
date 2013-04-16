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

#Function to find edge corresponding to matrix coordinate. t=tuple. l=level(index of m)
def findedge(t,l):
    


    

for i in range(ll):
    print m[i]

#Initial probability
inprob = []

for i in m[3].flat:
    inprob.append(i/m[ll-1].sum())

#Sample list
s = [0]*ll
p = [0]*ll


#Randomly choose first element in m[3] according to initial probabilities
s[ll-1] = np.random.choice(m[ll-1].size, p=inprob)
p[ll-1] = inprob[s[ll-1]]

print s[ll-1], p[ll-1]

#Iterate in reverse order
for i in range(ll-1, 0, -1):
    print 'i'
    print i
    
    prob = []
    idx = np.unravel_index(s[i], (n[i],n[i]))
    print idx

    #Need to calculate edges corresponding to idx so we know edge data that these tuples correspond to. function?


    
    
    for j in m[i-1].flat:
        print j
        print emission(GT[i],( , ))
        print haptrans()
        #prob.append((emission(GT[i], )*haptrans()*j) / m[i].flat[s[i]])
    
    print prob

    #print np.random.choice(globals()['a'+str(i)].size, p=prob)

print s

