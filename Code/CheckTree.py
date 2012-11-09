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

TestName=["DAG", "APERIODIC", "ROOT NODE", "LEAF"]
CheckList=[]
    
#DAG TEST: Check that the graph is a directed acyclic graph (DAG)
CheckList.append(is_directed_acyclic_graph(G))
        
       
#APERIODIC TEST: Check that graph is aperiodic
CheckList.append(is_aperiodic(G))


#ROOT NODE TEST: Check that the graph only contains one root node
PDic = nx.dfs_predecessors(G)
for i in range level

        #Check that there is only one root node




        #Check the number of leaves = no of haplotypes



        #Check that the number of nodes on each path is the length of the haplotypes no of levels = no of alleles\

        



def CheckTree(iterable):
    for i in range(iterable):
        if iterable[i] == False:
            return element
    return True


if CheckTree(CheckList) == True:
    print "TREE HAS PASSED CHECK"

else:
    
    print "TREE HAS FAILED CHECK DUE TO " + TestName[CheckTree(CheckLists)] + "TEST"
    
