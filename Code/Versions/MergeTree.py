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

from CreateTree import *

def plot_tree(G, name='plot.png'):
      
    #Create undirected graph from directed graph for visualization purposes           
    N=G.to_undirected()
    pyl.rcParams['figure.figsize'] = 10, 10
    
    #Clear plot to begin again
    plt.clf()
    plt.axis('off')
        
    #Set vizualization layout
    prog ='dot'
    pos = nx.drawing.graphviz_layout(G, prog)
    
    #Mark edges with labels corresponding to the weight of the edge
    edge_labels=dict([((u,v,),d['weight']) for u,v,d in N.edges(data=True)])
   
    
    
    #Create two edges lists for the different alleles
    allele_1 = [(u,v) for (u,v,d) in N.edges(data=True) if d['allele'] == '1']
    allele_2 = [(u,v) for (u,v,d) in N.edges(data=True) if d['allele'] == '2']
    #both_alleles = [(u,v) for (u,v,d) in N.edges(data=True) if d['weight'] >0.5]

    
    nodesize=[]
    for i in G.nodes():
        totalin = math.fsum([j[2]['weight'] for j in M.in_edges(i, data=True)])
        totalout = math.fsum([j[2]['weight'] for j in M.out_edges(i, data=True)])
        if totalin == totalout:
            nodesize.append(math.fabs(totalin))        
        else:
            nodesize.append(math.fabs(totalin-totalout))

            

    #Draw out nodes
    nx.draw_networkx_nodes(G,pos,node_size=nodesize, node_color='k')
    #Draw edges from edge lists of different alleles
    nx.draw_networkx_edges(N, pos, edgelist=allele_1, node_size=100, width=3, splines=True)
    nx.draw_networkx_edges(N, pos, edgelist=allele_2, node_size=100, width=6, style='dotted', splines=True)
    #nx.draw(N, pos, node_size=100, node_color='w', edge_color=edge_colours, width=4, with_labels=False)
    #Draw edge labels
    nx.draw_networkx_edge_labels(G,pos,ax=None,edge_labels=edge_labels)
    #Draw node labels
    #nx.draw_networkx_labels(G,pos,ax=None)

    

    #os.remove('plot.png')
    #Show plot in window
    #plt.show()

    #Save plot
    plt.savefig(name)

    
#Finds numerator of transition probability for given node and given allele. Returns 0 if allele is not found.
def gettransition(n,allele):
    t = 0

    #Each outgoing edge of the given node is checked to have the corresponding allele.
    for i in M.out_edges(n, data=True):
        if i[2]['allele'] == allele:
            t = i[2]['weight']
    
    return t

#Performs node test between a pair of nodes and all alleles which both nodes contain. Returns True or False.
def similarityscore(p,totala,totalb,threshold):
        
    #List of alleles which are on outgoing edges of both nodes being tested.
    allelelist = list(set(list([i[2]['allele'] for i in M.out_edges(p[0], data=True)]) + list([i[2]['allele'] for i in M.out_edges(p[1], data=True)])))

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
    totala = math.fsum([i[2]['weight'] for i in M.out_edges(a, data=True)])
    totalb = math.fsum([i[2]['weight'] for i in M.out_edges(b, data=True)])
   
    if totala == 0 and totalb == 0:
        return False
    
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
            for i in M.out_edges(q[0][0], data=True):
                x = i[2]['allele']
                y = None

                for j in M.out_edges(q[0][1], data=True):
                    
                    #If matching allele is found in outgoing edges of second node then for loop stops.
                    if j[2]['allele'] == x:
                        y = j
                        break
                                                
                    else:
                        continue
                    
                #Corresponding pairs of nodes with matching alleles are added to the queue.
                if y is not None:
                    q.append([i[1],y[1]])

            #Completed pair of nodes is removed from list          
            q.pop(0)
            
    #If all possible pairs are compared and are True, function returns True.
    return maxscore

#Function which returns the level of the given node
def whichlevel(n):    
    for i in range(len(mlevel)):
       if mlevel[i] <= n <= mlevel[i+1]:
	    return i+1
        
#Function which removes node from nlevel
def remove_node_from_nlevel(m):    
    for a in range(m,(haplolength++1)):
        nlevel[a] = nlevel[a]-1

#Inner merge function
def merge_inner(a,b):
	
    #List of edges directed out of nodes a and b.
    a_out = M.out_edges(a, data=True)
    b_out = M.out_edges(b, data=True)

    #Iterate through outgoing edges of node a
    for i in a_out:
        x = i[2]['allele']
        y = None
        
	#If outgoing edge of a has same allele as outgoing edge of b the index is set to y.
        for iterator in b_out:
            if iterator[2]['allele'] == x:
                y = iterator

	#If y exists, it is removed from list of outgoing edges from b
        if y is not None:
            b_out.remove(y)         
                       
            #get edge from b corresponding to i
            j = y
	    
            #Weight of edge in subtree of b is added to corresponding edge of a
            M.edge[i[0]][i[1]][0]['weight'] = M.edge[i[0]][i[1]][0]['weight'] + M.edge[j[0]][j[1]][0]['weight']
            
            #Recurse using child node
            merge_inner(i[1],j[1])

            #Remove edge and leaf node that is added to subgraph of a.
            M.remove_edge(y[0],y[1])
            M.remove_node(y[1])
            l=whichlevel(y[1])
            remove_node_from_nlevel(l) 

    #All outgoing edges of b which match one from a have been removing leaving only missing edges in list b_out
    for k in b_out:

	#Missing edges added to subgraph of a and remove the missing edge from graph
        M.add_weighted_edges_from([(a, k[1], k[2]['weight'])], allele = k[2]['allele'])
        M.remove_edge(k[0],k[1])
   
    
def mergenodes(a, b):
    #Call merge_inner
    merge_inner(a,b)

    #List of edges going into a
    b_in = M.in_edges(b, data=True)

    #Move all incoming edges of b in a. IMPROVEMENT = CAN REPLACE EDGE RATHER THAN ADD THEN DELETE
    for i in b_in:
        M.add_weighted_edges_from([(i[0], a, i[2]['weight'])], allele = i[2]['allele'])
        M.remove_edge(i[0],i[1])

    #Remove node b	
    M.remove_node(b)
    l=whichlevel(b)
    remove_node_from_nlevel(l)

   
               
    #Relabel node names which are above the level of the node being added.
    mapping=dict(zip(M.nodes(),range(1,nlevel[haplolength]+1)))
    nx.relabel_nodes(M,mapping,copy=False)

      
#Set mlevel to equal level
mlevel = list(level)
nlevel = list(mlevel)


#Create copy of G into multidirected graph
M = nx.MultiDiGraph()
M.add_nodes_from(G)
M.add_edges_from(G.edges_iter(data=True))
plot_tree(M, name='before_merge.png')

 
#Iterate through levels of M starting at second level.
for i in range(1, haplolength+1):
  
    nn = mlevel[i] - mlevel[i-1]
    
    
    #If there is only 1 node in graph G level i then carry on to next level. 
    if nn == 1:
        print 'option1'
        continue

    #If there is only 2 nodes in level. Test similarity between the pair of nodes.
    elif nn == 2:

              
        if mergetest(mlevel[i],mlevel[i]-1) == False:
            continue
        
        else:
            mergenodes(mlevel[i],mlevel[i]-1)
            mlevel = nlevel
            

    #If there is more than one node in level, test each pair of nodes in each level and merge the lowest scoring
    #Could make this more efficient by storing scores between nodes for each level so that after first node, merges are not repeated. ie. test between 2 nodes that you know is false.
    else:

        merge = True
        while merge == True:
            merge = False
            levelmin = 1000
            minj = 0
            mink = 0

            #Calculate merge score for each pair of nodes
            for j in range(mlevel[i-1]+1, mlevel[i]+1):
                for k in range(j, mlevel[i]+1):
                    
                    
                    if j != k:
                        testscore = mergetest(j,k)
                        
                       
                        if testscore == False:
                            continue
                        else:
                            if testscore < levelmin:
                                levelmin = testscore
                                minj = j
                                mink = k
                                            
            if levelmin != 1000:
                mergenodes(minj,mink)
                mlevel = nlevel
                merge = True


#Merge all nodes in final level
#Iterate through all outgoing edges of all the nodes on penultimate level
#Move edges to go to the first node of the final level
#Delete all other nodes after the first node of the final level.



endnode = mlevel[haplolength-1] + 1
for i in range(mlevel[haplolength-2]+1,mlevel[haplolength-1]+1):
    for j in M.out_edges(i, data=True):
        if j[1] != endnode:
            M.add_weighted_edges_from([(j[0], endnode, j[2]['weight'])], allele = j[2]['allele'])
            M.remove_edge(j[0],j[1])

for i in range(mlevel[haplolength-1]+2, mlevel[haplolength]+1):
    M.remove_node(i)
    mlevel[haplolength]=endnode

            


        
plot_tree(M, name='after_merge.png')

