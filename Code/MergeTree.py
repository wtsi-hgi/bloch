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
import pygraphviz as pgv
import operator
import os

from CheckInput import *

#Set mlevel to equal level
mlevel = level

#Create copy of G into multidirected graph
M = nx.MultiDiGraph()

M.add_nodes_from(G)
M.add_edges_from(G)

#Iterate through levels of M starting at second level.
for i in range(1, haplolength+1):

    nn = mlevel[i] - mlevel[i-1]
    
    #If there is only 1 node in graph G level i then carry on to next level. 
    if nn == 1:
        continue

    else if nn == 2:
        if nodepairtest(p) = True:
            

    #If there is more than one node in level, carry out merge test 
    else:
        
        
        
        
