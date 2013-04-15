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

  

print dir()
print a1
print a2
print a3
print a4

print a1.shape[0]

print a1.flatten()

inprob = []

for i in a4.flat:
    prob.append(i/a4.sum())

print inprob

#print numpy.version.version
#print np.random.choice(n4*n4, p=inprob)

s4 = 3

print np.unravel_index(choice, (n4,n4))

for i in range(ll-1, 0, -1):
    print 'i'
    print i
    prob = []
    
    for j in globals()['a'+str(i)].flat:

        prob.append(emission(GT[i-1],allele)*haptrans(e,d)*j/globals()['a'+str(i+1)].flat[globals()['s'+str(i+1)]])
    
    print prob

    #print np.random.choice(globals()['a'+str(i)].size, p=inprob)


#Collapse matrices into 1D array numpy.flatten function put.

#Initiation: Choose elements

#numpy.random.choice
