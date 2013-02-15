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
import math

#from CheckInput import *

#Splitting function takes a node as a variable and splits node according to haplotypes attached to the node.
def nodesplit(G,n,l,level):
    
    #For each haplotype   
    for i in range(len(G.node[n]['haplotype'])):
        #If edge exists, append frequency and haplotype suffix to the lists on next node.
        if G.node[n]['haplotype'][i][0] in [(edata['allele']) for u,v,edata in G.out_edges(n,data=True) if 'allele' in edata]:
            G.node[v]['haplotype'].append(G.node[n]['haplotype'][i][1:])
            G.node[v]['frequency'].append(G.node[n]['frequency'][i])

        #If edge does not exist, add edge.
        else:
            G.add_edge(n,level[l+1]+1,allele=G.node[n]['haplotype'][i][0])
            G.add_node(level[l+1]+1,haplotype=[G.node[n]['haplotype'][i][1:]],frequency=[G.node[n]['frequency'][i]])

            for a in range(l+1,(len(G.node[n]['haplotype'][0])+1)+l):
                level[a] = level[a] + 1      

#Function which tests if two nodes are similar enough to merge. Returns True or False.
def mergetest(a,b):
    print "a,b"
    print a,b
    #Score threshold and haplotype count for nodes a and b are calculated.
    totala = sum(G.node[a]['frequency'])
    totalb = sum(G.node[b]['frequency'])
    threshold = math.sqrt(math.pow(totala, -1) + math.pow(totalb, -1))
    maxscore = 0
    alevel = [1]*(len(G.node[a]['haplotype'][0])+1)
    blevel = [1]*(len(G.node[a]['haplotype'][0])+1)
    
    A=nx.DiGraph()
    A.clear()
    A.add_node(1,haplotype=G.node[a]['haplotype'],frequency=G.node[a]['frequency'])
    B=nx.DiGraph()
    B.clear()
    B.add_node(1,haplotype=G.node[b]['haplotype'],frequency=G.node[b]['frequency'])
    #for i in range(len(G.node[a]['haplotype'][0])-1 ):
    print 'Tree A nodes'
    for i in A.nodes_iter(data=True):
        print i
    print 'alevel'
    print alevel
    print 'Tree B nodes'
    for i in B.nodes_iter(data=True):
        print i
    print 'blevel'
    print blevel
        
    q = [[1,1,0]]

    while len(q) != 0:
        print "nodesplit A and B"
        print q[0][0],q[0][1],q[0][2]
        nodesplit(A,q[0][0],q[0][2],alevel)
        nodesplit(B,q[0][1],q[0][2],blevel)

        print 'Tree A edges'
        for i in A.edges(data=True):
            print i
        print 'Tree A nodes'
        for i in A.nodes_iter(data=True):
            print i
        print 'alevel'
        print alevel

        print 'Tree B edges'
        for i in B.edges(data=True):
            print i
        print 'Tree B nodes'
        for i in B.nodes_iter(data=True):
            print i
        print 'blevel'
        print blevel
        
        dicta = dict((edata['allele'],sum(A.node[v]['frequency'])) for u,v,edata in A.out_edges(q[0][0], data=True))
        dictb = dict((edata['allele'],sum(B.node[v]['frequency'])) for u,v,edata in B.out_edges(q[0][0], data=True))
        
        seta = set(dicta.keys())
        setb = set(dictb.keys())
        print "seta, setb"
        print seta, setb
        
        for i in list(seta.symmetric_difference(setb)):
            if i in dicta:
                print "in dicta"
                score = math.fabs(float(dicta[i])/float(totala))
                if score > threshold:
                    print "merge test false: allele not in both outgoing edges"
                    return False
                elif score > maxscore:
                    maxscore = score

            if i in dictb:
                print "in dictb"
                score = math.fabs(float(dictb[i])/float(totalb))
                if score > threshold:
                    print "merge test false: allele not in both outgoing edges"
                    return False
                elif score > maxscore:
                    maxscore = score
                   
        for i in list(seta.intersection(setb)):
            score = math.fabs(float(dicta[i])/totala - (float(dictb[i])/totalb))
            if score > threshold:
                print "merge test false: allele in both outgoing edges"
                return False
                
            else:
                print "edge test true"
                if score > maxscore:
                    maxscore = score
                print "i"
                print i
               
                print A.out_edges(q[0][0],data=True)
                a = [v for u,v,edata in A.out_edges(q[0][0],data=True) if edata['allele']==i]
                print a
                
                b = [v for u,v,edata in A.out_edges(q[0][0],data=True) if edata['allele']==i]
                print b

                    
                #for u,v,edata in B.out_edges(q[0][1],data=True) if edata['allele']=='i':
                #print u,v,edata
                if len(A.node[a[0]]['haplotype'][0]) != 0:
                    print "q.append"
                    print a[0],b[0],q[0][2]+1
                    q.append([a[0] ,b[0] ,q[0][2]+1])
                print "q is"
                print q
            
        q.pop(0)
        print "End: q is"
        print q
        
    print "return maxscore"  
    return maxscore

def mergenodes(a,b):
    print 'Edges in Graph'
    for i in G.edges(data=True):
        print i
    print 'Nodes in Graph'
    for i in G.nodes_iter(data=True):
        print i
    print 'glevel'
    print glevel

    for i in G.node[a]['haplotype']:
        if i in G.node[b]['haplotype']:
            print "yes"

    seta=[i for i in G.node[a]['haplotype']]
    setb=[i for i in G.node[b]['haplotype']]
    print seta
    print setb

    for i in list(seta.symmetric_difference(setb)):
        print i
        #only if summetric_difference exists!

    for i in list(seta.intersection(setb)):
        print i

    #add frequencies on to same haplotypes
    #move edge to second node to the first node.
    #delete last node and edit level.
   

def merge(l):
    nn = glevel[l] - glevel[l-1]
    if nn == 1:
        print "one node in level so no merge can be performed"
        pass
    
    #If there is only 2 nodes in level. Test similarity between the pair of nodes.
    elif nn == 2:
        print "2 nodes in level"
               
        if mergetest(glevel[l],glevel[l]-1) == False:
            pass
        
        else:
            print "nodes need to be merged"
            mergenodes(level[i],level[i]-1)

    else:
        print "more than two nodes in level"        
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
                        print "j,k,testscore"
                        print j,k,testscore
                 
                        if testscore == False:
                            continue
                        else:
                            if testscore < levelmin:
                                levelmin = testscore
                                minj = j
                                mink = k
                                            
            if levelmin != 1000:
                print "merge nodes"
                print minj,mink
                mergenodes(minj,mink)
                merge = True

haplotype = ['1111', '1112', '1122', '1221', '2111', '2112', '2122']
frequency = [21, 79, 95, 116, 25, 112, 152]
haplolength=4
haplonum=7

glevel=[1]*(haplolength+1)
G=nx.MultiDiGraph()

#Add start node 1
G.add_node(1,haplotype=haplotype,frequency=frequency)

#Split node 1 of graph G.
print "nodesplit(G,1,0,glevel)"
nodesplit(G,1,0,glevel)




#Carry out merge part on all nodes in level 1
print "merge(1)"
merge(1)

#for i in range(glevel[l-1]+1,glevel[l]+1):
#nodesplit(G,i,l)


print "nodesplit(G,2,i,glevel) i=[1,2]"
for i in range(glevel[0]+1,glevel[1]+1):
    nodesplit(G,i,1,glevel)

print 'Edges in Graph'
for i in G.edges(data=True):
    print i
print 'Nodes in Graph'
for i in G.nodes_iter(data=True):
    print i
print 'glevel'
print glevel


print "merge(2)"
#merge(2)
mergenodes(4,6)

