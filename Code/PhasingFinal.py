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

#Extract genotypes from data
with open('/Users/mp18/Documents/bloch/Data/Table_1', 'rb') as f:
    reader = csv.reader(f, delimiter='\t',skipinitialspace = True)
    #Create list of tuples which contain genotypes of each individual
    GT = zip(*reader)

#Remove first tuple of position names
del GT[0]

#Create input for tree algorithm
haplotypes = {}

#Randomly assign phase for each sample
for i in GT:
    a=''
    b=''
    for j in i:
        x=j[0]
        y=j[-1]
        rand = randrange(0,2)
        if rand == 0:
            a+=x
            b+=y
        else:
            a+=y
            b+=x
    
    if a in haplotypes:
        haplotypes[a] += 1
    else:
        haplotypes[a] = 1

    if b in haplotypes:
        haplotypes[b] += 1
    else:
        haplotypes[b] = 1


#Haplotype length
hlength = len(a)-1

def treealgorithm(h):      

    #Number of different haplotypes
    hnum = len(h)
    
    #Function which takes node and splits node according to haplotype dictionary attached to the node
    def nodesplit(n, l):
        #Loop over keys in node's haplotype dictionary
        for key, value in G.node[n]['hap'].iteritems():
        
            #If there exists an edge which corresponds to first character of key.
            if key[0] in [(edata['a']) for u,v,k,edata in G.out_edges(n,data=True, keys=True)]:
                #Add value to weight variable on edge
                G.edge[u][v][k]['weight'] += value
                #Add haplotype suffix and value to dictionary of connecting node
                G.node[v]['hap'].update({key[1:]:value})

            #If there does not exist an edge.
            else:
                m = gl[l+1]+1
                #Create edge and label it with first character of key
                G.add_edge(n,m,a=key[0],weight=value)
                #Assign dictionary of suffix of that key to node that has been created.
                G.add_node(m, hap={key[1:]:value})

                #Add node to gl
                for i in range(l+1, hlength+1):
                    gl[i] += 1


    #Merge function carries out pairwise test between all nodes on given level and merges the lowest scoring nodes.
    #Cycle is repeated until no more merges can be made on the give level.
    def merge(l):

        #Set variable nn for number of nodes
        nn = gl[l] - gl[l-1]
    
        #If only one node on level.
        if nn == 1:
            #No need to carry out merge function so returns.
            return
    
        #If there is only 2 nodes in level. 
        elif nn == 2:
            #No need to carry out merge if nodes fail mergetest
            if mergetest(gl[l],gl[l]-1) == False:
                return
            #Merge nodes
            else:
                mergenodes(gl[l],gl[l]-1)

        #If there are more than 2 nodes in a level, all pairs of nodes are compared and lowest scoring is merged.
        else:
            merge = True
            while merge == True:
                merge = False
                levelmin = -1
                minj = 0
                mink = 0
            
                #Calculate merge score for each pair of nodes
                for j in range(gl[l-1]+1, gl[l]+1):
                    for k in range(j, gl[l]+1):                    
                        if j != k:                        
                            testscore = mergetest(j,k)
                            if testscore == False:
                                continue
                            else:                            
                                #Details of lowest scoring pair are stored
                                if testscore < levelmin:
                                    levelmin = testscore
                                    minj = j
                                    mink = k
                                
                #The lowest scoring pair of nodes are merged
                if levelmin != -1:                
                    mergenodes(minj,mink)                
                    #Whenever a merge has taken place, merge is set to True so that the loop is repeated
                    merge = True






    
    #Create networkx MultiGraph
    G=nx.MultiDiGraph()

    #Add start node 1 and split the first node.
    G.add_node(1,hap=h)
    nodesplit(1,0)
    merge(1)
   

    
    return G

#Input into tree algorithm
G = treealgorithm(haplotypes)

#Carry out forward algorithm and backward sampling once conditional on each genotype. Save output of backward sampling

#Reverse marker order. Use as input to next iteration.


