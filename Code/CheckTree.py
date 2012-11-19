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

from CreateTree import *

TestName=[]
CheckList=[]
    
#DAG TEST: Check that the graph is a directed acyclic graph (DAG)
CheckList.append(is_directed_acyclic_graph(G))
TestName.append("DAG")
        
#APERIODIC TEST: Check that graph is aperiodic
CheckList.append(is_aperiodic(G))
TestName.append("APERIODIC")

#ROOT NODE TEST: Check that the first level only contains one root node
#ie. Level[1] begins with number 2
if level[1] == 3:
    CheckList.append(True)
else:
    CheckList.append(False)

TestName.append("ROOT NODE")

#LEAVES TEST: Check that the number of nodes on the last level is equal to the
#number of unique haplotypes.
hapnum = len(list(set(haplotype)))

if hapnum == (level[len(level)] - level[len(level)-1]):
    CheckList.append(True)
else:
    CheckList.append(False)

TestName.append("LEAVES")

#EDGES TEST: Check that each edge is connected between nodes which are on consecutive levels

def whichlevel(node):
    for i in range(len(level)):
        if level[i] < node <= level[i+1]:
            return 



for i in nx.bfs_edges(G, 1):
    




TestName.append("EDGES")







        #Check that the number of nodes on each path is the length of the haplotypes no of levels = no of alleles\

        
        #Check there exists a path between root node and all leaf nodes
        #Check length of all paths is the same + same as length of haplotype

def CheckTree(CheckList):
    for i in range(CheckList):
        if CheckList[i] == False:
            return element
    return True


if CheckTree(CheckList) == True:
    print "TREE HAS PASSED CHECK"

else:
    
    print "TREE HAS FAILED CHECK DUE TO " + TestName[CheckTree(CheckList)] + "TEST"
    
