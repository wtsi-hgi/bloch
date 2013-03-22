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


from v3 import G, glevel
print "Nodes"
for i in G.nodes(data=True):
    print i
    
print "Edges"
for i in G.edges(data=True):
    print i

print "Adjacency List"
for i in G.adjacency_iter():
    print i

print G[6].keys()

print "glevel"
print glevel
ll = len(glevel) -1


GT = (('?','?'),(1,1),(1,2),(1,2))

#Initial state probabilities
def hapinitial(allele):
    for i in G.out_edges(1, data=True):
        
        if i[2]['allele'] == allele:
            return sum(G.node[i[1]]['frequency'])/sum(G.node[1]['frequency'])

def dipinitial(a,b):
    return hapinitial(a)*hapinitial(b)

#Emission state probabilities
def emission(gt,s):
    if gt.count('?') == 1:
        print '1'
    elif gt.count('?') == 2:
        print '2'
        

        
    else:
        if set(gt) == set(s):
            return 1
        else:
            return 0

#Initiation
for i in G[1].keys():
    for j in G[1][i].keys():
        print G[1][i][j]

for i in G.out_edges(1):
    for j in G.out_edges(1):
        print i,j

print dipinitial(1,2)
print emission(('?','?'),(1,1))

#Induction


