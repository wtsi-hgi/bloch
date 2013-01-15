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
import networkx as nx
import sys
import matplotlib.pyplot as plt
import pylab as pyl
import pygraphviz as pgv
import operator
import os
import math
from networkx.utils import is_string_like
from collections import deque


def edgetest(a,b,allele):
    totala = math.fsum([i[2]['weight'] for i in A.out_edges(a, data=True)])
    totalb = math.fsum([i[2]['weight'] for i in A.out_edges(b, data=True)])
    threshold = math.sqrt(math.pow(totala, -1) + math.pow(totalb, -1))

    print 'threshold'
    print threshold

    for i in A.out_edges(a, data=True):
        if i[2]['allele'] == allele:
            aw = i[2]['weight']
        else:
            aw = 0

    for i in A.out_edges(b, data=True):
        if i[2]['allele'] == allele:
            bw = i[2]['weight']
        else:
            bw = 0

    score = math.fabs(aw/totala - bw/totalb)

    print 'score'
    print score

    if score > threshold:
        return False
    else:
        return True

def mergetest(a,b):
    
    allelelist = []
    
    for i in A.out_edges(2, data=True):
        if i[2]['allele'] not in allelelist:
            allelelist.append(i[2]['allele'])

    test = True
    for i in allelelist:
        
        test = edgetest(2,3,i)
    
        if test is False:
            return False





    #for i in A.get_edge_data(a,b):
        

    #qa = deque()
    #qb = deque()

    #qa.append(a)
    #qb.append(b)

    #while len(qa)]]['

#Create test graph manually
A = nx.MultiDiGraph()
A.add_weighted_edges_from([(1,2,311),(2,4,195),(4,7,100),(7,12,21),(3,6,289),(6,10,137),(10,16,25)], allele=1)
A.add_weighted_edges_from([(7,13,79),(4,8,95),(8,14,95),(2,5,116),(5,9,116),(9,15,116),(1,3,289),(10,17,112),(6,11,152),(11,18,152)], allele=2)


#Set mlevel to
mlevel = [1,3,6,11,18]
haplolength = 4

print mergetest(2,3)







