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

#Create list of haplotypes
haplotype = []

#Create list of haplotype frequencies
frequency = []

#Open data file in read mode
with open('/Users/mp18/Documents/bloch/Data/Table_3', 'r') as f:
    data = f.read()
    
#Close file
f.close()

#Seperate data file line and add each item to corresponding list.
count = 0

for word in data.split():    
    if count==0:
        #Haplotypes stored as strings and appended to the list haplotype.
        haplotype.append(word)
        count = 1
        continue
    
    elif count==1:        
        try:
            #Frequencies stored as integers and appended to the list frequency
            num = int(word)
            
        except ValueError:
            print "Error: Frequency data contains a non-integer"
        
        frequency.append(num)
        count = 0
        continue

#Set length of first haplotype to equal haplolength
haplolength = len(haplotype[0])

#Set number of haplotypes to equal haplonum
haplonum = len(haplotype)


#For each haplotype
for i in range(haplonum):

    #Checks that each haplotype in each iteration is the same length. If not, error is printed.
    if haplolength != len(haplotype[i]):
        print ("Error: Haplotype number " + str(i) + " is not the same length") 
        break


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
    edge_labels={}
    for i in G.edges_iter(data=True, keys=True):
        
        if len(i[3]) == 1:
            edge_labels[(i[0],i[1])] = sum(G.node[i[1]]['frequency'])
        else:
            if i[2] == 0:
                edge_labels[(i[0],i[1])] = i[3]['weight']
               
            else:
                weight = edge_labels[(i[0],i[1])] + i[3]['weight']  
                edge_labels[(i[0],i[1])] = weight  
    
    #Create two edges lists for the different alleles
    allele_1 = [(u,v) for (u,v,d) in N.edges(data=True) if d['allele'] == 1]
    allele_2 = [(u,v) for (u,v,d) in N.edges(data=True) if d['allele'] == 2]
    allele_3 = [(u,v) for (u,v,d) in N.edges(data=True) if d['allele'] == 3]
        
    nodesize=[]
    for i in G.nodes():
        nodesize.append(sum(G.node[i]['frequency']))

    #Draw out nodes
    nx.draw_networkx_nodes(G,pos,node_size=nodesize, node_color='k')
    #Draw edges from edge lists of different alleles
    nx.draw_networkx_edges(N, pos, edgelist=allele_1, node_size=100, width=3, splines=True)
    nx.draw_networkx_edges(N, pos, edgelist=allele_2, node_size=100, width=6, style='dotted', splines=True)
    nx.draw_networkx_edges(N, pos, edgelist=allele_3, node_size=100, width=6, style='dashed', splines=True)

    #Draw edge labels
    nx.draw_networkx_edge_labels(G,pos,ax=None,edge_labels=edge_labels)
    #Draw node labels
    #nx.draw_networkx_labels(G,pos,ax=None)    

    #Save plot
    plt.savefig(name)

#Splitting function takes a node as a variable and splits node according to haplotypes attached to the node.
def nodesplit(G,n,l,level):
    
    #For each haplotype   
    for i in range(len(G.node[n]['haplotype'])):
        #If edge exists, append frequency and haplotype suffix to the lists on next node.
        if int(G.node[n]['haplotype'][i][0]) in [(edata['allele']) for u,v,edata in G.out_edges(n,data=True) if 'allele' in edata]:
            G.node[v]['haplotype'].append(G.node[n]['haplotype'][i][1:])
            G.node[v]['frequency'].append(G.node[n]['frequency'][i])

        #If edge does not exist, add edge.
        else:
            G.add_edge(n,level[l+1]+1,allele=int(G.node[n]['haplotype'][i][0]))
            G.add_node(level[l+1]+1,haplotype=[G.node[n]['haplotype'][i][1:]],frequency=[G.node[n]['frequency'][i]])

            for a in range(l+1,(len(G.node[n]['haplotype'][0])+1)+l):
                level[a] = level[a] + 1      

#Function which tests if two nodes are similar enough to merge. Returns True or False.
def mergetest(a,b):  
    
    #Formula components are calculated
    totala = sum(G.node[a]['frequency'])
    totalb = sum(G.node[b]['frequency'])
    threshold = math.sqrt(math.pow(totala, -1) + math.pow(totalb, -1))
    maxscore = 0
    alevel = [1]*(len(G.node[a]['haplotype'][0])+1)
    blevel = [1]*(len(G.node[a]['haplotype'][0])+1)    

    #Subgraphs A and B are created in order to calculate similarity score
    A=nx.DiGraph()
    A.clear()
    A.add_node(1,haplotype=G.node[a]['haplotype'],frequency=G.node[a]['frequency'])
    B=nx.DiGraph()
    B.clear()
    B.add_node(1,haplotype=G.node[b]['haplotype'],frequency=G.node[b]['frequency'])
        
    q = [[1,1,0]]

    while len(q) != 0:
        
        #Splits nodes subgraphs
        nodesplit(A,q[0][0],q[0][2],alevel)
        nodesplit(B,q[0][1],q[0][2],blevel)
        
        dicta = dict((edata['allele'],sum(A.node[v]['frequency'])) for u,v,edata in A.out_edges(q[0][0], data=True))
        dictb = dict((edata['allele'],sum(B.node[v]['frequency'])) for u,v,edata in B.out_edges(q[0][0], data=True))
        
        seta = set(dicta.keys())
        setb = set(dictb.keys())

        #Iterates through alelles which are not a member of both node's outgoing edges and tests score against threshold
        for i in list(seta.symmetric_difference(setb)):
            if i in dicta:
                score = math.fabs(float(dicta[i])/float(totala))
                if score > threshold:
                    return False
                elif score > maxscore:
                    maxscore = score

            if i in dictb:
                score = math.fabs(float(dictb[i])/float(totalb))
                if score > threshold:    
                    return False
                elif score > maxscore:
                    maxscore = score

        #Iterates through alleles which are a member of both node's outgoing edges
        for i in list(seta.intersection(setb)):
            score = math.fabs(float(dicta[i])/totala - (float(dictb[i])/totalb))
            if score > threshold:   
                return False                
            else:
                if score > maxscore:
                    maxscore = score               
                a = [v for u,v,edata in A.out_edges(q[0][0],data=True) if edata['allele']==i]                               
                b = [v for u,v,edata in B.out_edges(q[0][0],data=True) if edata['allele']==i]
                
                #Adds next nodes to list if edges has passed edge test            
                if len(A.node[a[0]]['haplotype'][0]) != 0:
                    q.append([a[0],b[0],q[0][2]+1])
                    
        q.pop(0)
  
    return maxscore

#Function merges 2 nodes
def mergenodes(a,b):
    global G

    #Add weight attribute to merged edge so that information is not lost
    for i in G.in_edges(a, data=True, keys=True):
        G.add_edge(i[0],i[1],key=i[2],allele=i[3]['allele'],weight=sum(G.node[a]['frequency']))
            
    #Iterates through list of haplotypes on second node
    for i in range(len(G.node[b]['haplotype'])):
        #Add weights of haplotypes of b to a       
        if G.node[b]['haplotype'][i] in [j for j in G.node[a]['haplotype']]:
            index = [j for j in G.node[a]['haplotype']].index(G.node[b]['haplotype'][i])
            G.node[a]['frequency'][index] = G.node[a]['frequency'][index] + G.node[b]['frequency'][i]
        else:
            G.node[a]['haplotype'].append(G.node[b]['haplotype'][i])
            G.node[a]['frequency'].append(G.node[b]['frequency'][i])

    #Move all incoming edges of b to a
    for i in G.in_edges(b, data=True):
        G.add_edge(i[0],a,allele=i[2]['allele'],weight=sum(G.node[b]['frequency']))
        G.remove_edge(i[0],i[1])   

    G.remove_node(b)
    
    #Node deleted on glevel
    for i in glevel:
        if b <= i:
            glevel[glevel.index(i)] = glevel[glevel.index(i)] - 1

    #Relablel nodes so that they are consecutive integers
    G = nx.convert_node_labels_to_integers(G, first_label=1, ordering="sorted")

#Merge function carries out pairwise test between all nodes on each level and merges the lowest scoring
#nodes. This repeats in a cycle on each level until no more merges can be made
def merge(l):
    
    nn = glevel[l] - glevel[l-1]
    
    #If only one node on level, there is no need to carry out merge function.
    if nn == 1:     
        return
    
    #If there is only 2 nodes in level. Test similarity between the pair of nodes.
    elif nn == 2:               
        if mergetest(glevel[l],glevel[l]-1) == False:
            pass        
        else:
            mergenodes(glevel[i],glevel[i]-1)

    #If there are more than 2 nodes in a level, all pairs of nodes are compared and lowest scoring is merged.
    else:
        merge = True
        while merge == True:
            merge = False
            levelmin = 1000
            minj = 0
            mink = 0
            #Calculate merge score for each pair of nodes
            for j in range(glevel[l-1]+1, glevel[l]+1):
                for k in range(j, glevel[l]+1):                    
                    if j != k:
                        
                        testscore = mergetest(j,k)
                        if testscore == False:
                            continue
                        else:
                            #Details of lowest scoring pair are kept.
                            if testscore < levelmin:
                                levelmin = testscore
                                minj = j
                                mink = k
                                
            #All pairs of nodes on current level are retested if a merge has taken place.
            if levelmin != 1000:
                mergenodes(minj,mink)
                merge = True

#Create list of nodes in each level for graph G. Each position is a cumulative count
#of the number of nodes on the level and the previous levels.
glevel=[1]*(haplolength+1)

#Create networkx MultiGraph
G=nx.MultiDiGraph()

#Add start node 1 and split the first node.
G.add_node(1,haplotype=haplotype,frequency=frequency)
nodesplit(G,1,0,glevel)

#Iterate through each successive level.
for i in range(1,haplolength):

    #Carry out the merge function for each level.
    merge(i)
    

    #Split each node on each level
    for j in range(glevel[i-1]+1,glevel[i]+1):
        nodesplit(G,j,i,glevel)

#Set endnode as first node in last level
endnode = glevel[haplolength-1] + 1

#Add weight attribute to all incoming edges so that information is not lost
for i in G.in_edges(endnode, data=True, keys=True):
    G.add_edge(i[0],i[1],key=i[2],allele=i[3]['allele'],weight=sum(G.node[i[1]]['frequency']))

#Add edges ending in last level are moved to connect to endnode
for i in range(glevel[haplolength-2]+1,glevel[haplolength-1]+1):   
    for j in G.out_edges(i, data=True):        
        if j[1] != endnode:
            G.add_edge(j[0],endnode,allele=j[2]['allele'],weight=sum(G.node[j[1]]['frequency']))
            G.node[endnode]['frequency'][0]=G.node[endnode]['frequency'][0]+G.node[j[1]]['frequency'][0]
            G.remove_edge(j[0],j[1])

#Remove all nodes which are labelled above the endnode
for i in range(glevel[haplolength-1]+2, glevel[haplolength]+1):
    G.remove_node(i)

#Set glevel so it is correct    
glevel[haplolength]=endnode
    
#Plot resulting tree
plot_tree(G, name='v3.png')

