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

#Finds numerator of transition probability for given node and given allele. Returns 0 if allele is not found.
def gettransition(n,allele):
    t = 0

    #Each outgoing edge of the given node is checked to have the corresponding allele.
    for i in A.out_edges(n, data=True):
        if i[2]['allele'] == allele:
            t = i[2]['weight']
    
    return t

#Performs node test between a pair of nodes and all alleles which both nodes contain. Returns True or False.
def similarityscore(p,totala,totalb,threshold):
        
    #List of alleles which are on outgoing edges of both nodes being tested.
    allelelist = list(set(list([i[2]['allele'] for i in A.out_edges(p[0], data=True)]) + list([i[2]['allele'] for i in A.out_edges(p[1], data=True)])))

    #Score is calculated using each allele in list
    for i in allelelist:
        score = math.fabs(gettransition(p[0],i)/totala - gettransition(p[1],i)/totalb)

        #Function returns True if score is below the threshold.
        if score > threshold:
            return False
        else:
            return score

#Function which tests if two nodes are similar enough to merge. Returns True or False.
def mergetest(a,b):

    #Score threshold and haplotype count for nodes a and b are calculated.
    totala = math.fsum([i[2]['weight'] for i in A.out_edges(a, data=True)])
    totalb = math.fsum([i[2]['weight'] for i in A.out_edges(b, data=True)])    
    threshold = math.sqrt(math.pow(totala, -1) + math.pow(totalb, -1))
    maxscore = 0

    #Queue list is created for breadth-first search.
    q = [[a,b]]

    #Similarity score is calculated for the first item of the queue as long as there are elements in the queue.
    while len(q) != 0:

        score = similarityscore(q[0],totala,totalb,threshold)
        
        #If a pair of nodes are not similar enough mergetest returns False immediately.
        if score is False:
            return False

        else:
            if score > maxscore:
                maxscore = score
            
            #Each allele of outgoing edges of first node is compared to alleles from outgoing edges of second node.
            for i in A.out_edges(q[0][0], data=True):
                x = i[2]['allele']
                y = None

                for j in A.out_edges(q[0][1], data=True):
                    
                    #If matching allele is found in outgoing edges of second node then for loop stops.
                    if j[2]['allele'] == x:
                        y = j
                        break
                                                
                    else:
                        continue
                    
                #Corresponding pairs of nodes iwth matching alleles are added to the queue.
                if y is not None:
                    q.append([i[1],y[1]])

            #Completed pair of nodes is removed from list          
            q.pop(0)
            
    #If all possible pairs are compared and are True, function returns True.
    return maxscore
       

#Create test graph manually
A = nx.MultiDiGraph()
A.add_weighted_edges_from([(1,2,311),(2,4,195),(4,7,100),(7,12,21),(3,6,289),(6,10,137),(10,16,25)], allele=1)
A.add_weighted_edges_from([(7,13,79),(4,8,95),(8,14,95),(2,5,116),(5,9,116),(9,15,116),(1,3,289),(10,17,112),(6,11,152),(11,18,152)], allele=2)

#Set mlevel to
mlevel = [1,3,6,11,18]
haplolength = 4

print mergetest(4,6)



