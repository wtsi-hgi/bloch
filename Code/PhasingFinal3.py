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
import csv
import networkx as nx
import random
import numpy as np
import itertools
import math
import sys
import pickle
import cPickle
import string
import functools
import argparse
import time

#Tree class. T = Tree. tnodes = list of nodes. d = direction. 0 = forward. 1 = backwards. c = node count.
class Tree:
    def __init__(self, direction):
        #Initialize variables
        self.T = nx.MultiDiGraph()
        self.tnodes = []
        self.d = direction
        self.c = 1        
        
        #Add first node
        self.T.add_node(1,level=0,weight=0)
        for i in range(hlength+2):
            self.tnodes.append([])
        self.tnodes[0].append(1)
        
    #Method which adds a genotype to the graph. Input of first iteration has already been randomly phased and imputed.
    #gt is a list of lists. Each marker corresponds to a list with 2 elements in.
    def add_genotype(self, gt):
        
        #Add single marker position to the graph. cnode = current node. a = allele. l = level
        def add_marker(cnode,a,l):
            #If current node has no outgoing edges
            if self.T.out_degree(cnode) == 0:
                #Increase node name counter by one
                self.c += 1
                
                #Add edge between current node and new node.
                self.T.add_edge(cnode,c,a=a,weight=1)
                self.T.add_node(c, level=l+1,weight=1)
                self.tnodes[l+1].append(c)
                
                #Return current node position
                return self.c
            
            else:
                #Search through all outgoing nodes to find corresponding allele
                for u,v,k,edata in self.T.out_edges(cnode,data=True,keys=True):     
                    #If corresponding allele found. Add weight to edge and node count.           
                    if edata['a'] == a:
                        self.T.edge[u][v][k]['weight'] += 1
                        self.T.node[v]['weight'] +=1
                        
                        #Return current node position
                        return v
                
                #If no match is found add new edge and node
                self.c +=1
                self.T.add_edge(cnode,self.c,a=m,weight=1)
                self.T.add_node(self.c,level=l+1,weight=1)
                self.tnodes[l+1].append(self.c)
                
                #Return current node position
                return self.c               
            
        #Weight added to first node
        self.T.node[1]['weight'] += 2   
             
        #Current nodes set for haplotypes a and b
        anode = 1
        bnode = 1
        
        #Iterate through each marker in the position and add to the graph
        for n, i in enumerate(gt):        
            anode = add_marker(anode, i[0], n)           
            bnode = add_marker(bnode, i[1], n)
        
    #Merge nodes on each level
    def merge(self):
        pass
        
        


#Beginning of code.  Print start time.
sys.stdout.write("Start\t"+str(time.clock())+"\n")
T = Tree(0)
