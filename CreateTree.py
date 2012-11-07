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

import networkx as nx
import sys
import matplotlib.pyplot as plt

#Create list of haplotypes
haplotype = []

#Create list of haplotype frequencies
frequency = []

#Open data file in read mode
with open('/Users/mp18/Documents/Python/Data/Table_3', 'r') as f:
    data = f.read()
f.close()

#Seperate data file line and add each item to corresponding list.
count = 0
for word in data.split():
    if count==0:
        #Haplotypes stored as strings.
        haplotype.append(word)
        count = 1
        continue
    elif count==1:
        #Frequencies stored as integers.
        num = int(word)
        frequency.append(num)
        count = 0
        continue

#Set length of first haplotype to equal haplolength
haplolength = len(haplotype[1])


#Set number of haplotypes to equal haplonum
haplonum = len(haplotype)


#Create list of levels. Each position is a cumulative count of the number
#of nodes on the level and the previous levels.
level = [0]*(haplolength+1)

#Set current node pointer
cnode = 0


#Create empty directed graph
G = nx.DiGraph()


#Function to add node to level count. l=level
def add_node_to_level(m):

    #If node is added to last level
    if m == haplolength:
        level[m] = level[m] + 1

    #If node is added to a level which is not the last
    else:
        for a in range(m,(haplolength+1)):
            level[a] = level[a] + 1
      
        
#Function to add node to current node c=cnode, h=haplotype, l=level
def add_edge_from_cnode(G,c,h,l):

    #If node is being created at last level
    if l == haplolength - 1:
        G.add_edge(c, (level[l+1]+1), weight=frequency[h], allele=haplotype[h][l])

    #If node is being created that is not at the last level
    else:
       
        #Relabel nodes in range(level[j+1]+1,level[j+2]+1) to one number above
        def mapping(x):
            if x > level[l+1]:
                return x+1
            else:
                return x

        nx.relabel_nodes(G,mapping,copy=False)
        
    
        #Add edge from cnode to level[j+1]+1
        G.add_edge(c, (level[l+1]+1), weight=frequency[h], allele=haplotype[h][l])
       
        print 'nodes'
        print G.edges(data=True)
    
    

#Function which checks if edge from node exists. n=node to be checked. a=allele.
def check_for_existing_edge(n, a):
    for i in G.successors_iter(n):
        if G[n][i]['allele'] == a:
            return i
    return 0




#For each haplotype
for i in range(haplonum):

    #Checks that each haplotype in each iteration is the same length
    if haplolength != len(haplotype[i]):
        print 'Error: Not all haplotypes are the same length'
        break
    
    #If all haplotypes are the same length
    else:
        
        cnode = 1

        #For each allele in haplotype
        for j in range(haplolength):
            
            #If no nodes have been created, add first edge between node 1 and node 2
            if level[haplolength] == 0:
                G.add_edge(1, 2, weight=frequency[i], allele=haplotype[i][j])
                add_node_to_level(0)
                add_node_to_level(1)
                cnode = 2

            #If current node has no outgoing edges.
            elif len(G.successors(cnode)) == 0:            
                add_edge_from_cnode(G,cnode,i, j)
                cnode = level[j+1]+1
                add_node_to_level(j+1)
                 
                                
                
            #If current node has at least one outgoing edge 
            else:
                snode = check_for_existing_edge(cnode, haplotype[i][j])
                
                
                #If no edges have the same allele
                if snode == 0:           
                    add_edge_from_cnode(G,cnode,i, j)
                    add_node_to_level(j+1)
                #If there exists an edge with the same allele
                else:

                    G[cnode][snode]['weight'] = G[cnode][snode]['weight'] + frequency[i]
                    cnode = snode

           
       
                    
                    
                
labels=nx.draw_networkx_labels(G,pos=nx.spring_layout(G))
edge_labels=nx.draw_networkx_edge_labels(G,pos=nx.spring_layout(G))
nx.draw_networkx(G)
plt.show()


           
	
	





