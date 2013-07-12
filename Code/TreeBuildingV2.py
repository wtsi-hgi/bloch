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

import argparse
import csv
import networkx as nx
import random
import sys
import time
import itertools
import numpy as np


def random_weighted_choice(weights):
    rnd = random.random() * sum(weights)
    for i, w in enumerate(weights):
        rnd -= w
        if rnd < 0:
            return i

def treebuild(GT):

    #Add marker return new marker
    def add_marker(c, m, l):
        if G.out_degree(c) == 0:
            n = G.order()+1
            G.add_edge(c,n,a=m,weight=1)
            G.add_node(n, level=l+1)
            gnodes[l+1].append(n)
            return n
        else:
            for u,v,k,edata in G.out_edges(c,data=True, keys=True):
                
                if edata['a'] == m:
                    G.edge[u][v][k]['weight'] += 1
                    return v

            n = G.order()+1
            G.add_edge(c,n,a=m,weight=1)
            G.add_node(n, level=l+1)
            gnodes[l+1].append(n)
            return n      
    
    def add_genotype(gt):
        #Set node pointers for adding both haplotypes
        anode = 1
        bnode = 1
        
        for n, i in enumerate(gt):
            if i[1] == "|":
                a = i[0]
                b = i[2]                
            elif i[1] == "/":
                rand = random.randrange(0,2)
                if rand == 0:
                    a=i[0]
                    b=i[2]
                else:
                    a=i[2]
                    b=i[0]         
            else:
                print i[1]
                print "Unknown symbol"

            if a == '.':
                a = alleles[n][random_weighted_choice(allelefreq[n])]
            if b == '.':
                b = alleles[n][random_weighted_choice(allelefreq[n])]          

            anode = add_marker(anode, a, n)           
            bnode = add_marker(bnode, b, n)

    def mergetest(a,b):
        pass

    def mergenodes(a,b):
        #Move all incoming edges of b to a
        for i in G.in_edges(b, data=True, keys=True):
            G.add_edge(i[0],a,a=i[3]['a'],weight=i[3]['weight'])
            G.remove_edge(i[0],i[1],key=i[2])

        q=[[a,b]]

        while len(q) != 0:

            da = [edata['a'] for u,v,edata in G.out_edges(a, data=True)]
            db = [edata['a'] for u,v,edata in G.out_edges(b, data=True)]     
            
            sa = set(da)           
            sb = set(db)

            print sa
            print sb
            
            #Iterates through alelles which elements in sb but not in sa
            for i in sb.difference(sa):                
                if i in da:
                    pass
                    

                if i in db:
                    pass
            
            #Iterates through alleles which are a member of both node's outgoing edges
            for i in sa.intersectionn(sb):
                pass




            q.pop(0)

        


        G.remove_node(b)


        
    
           
    gnodes = []
    #Create list of nodes on each level.
    for i in range(hlength+2):
        gnodes.append([])    

    #Create networkx MultiGraph
    G=nx.MultiDiGraph()
    
    #Add start node 1.
    G.add_node(1,level=0)    
    gnodes[0].append(1)    
    
    for i in GT:
        add_genotype(i)

    for j in G.nodes(data=True):
        print j

    for j in G.edges(data=True, keys=True):
        print j

    for i in gnodes:
        if len(i) == 1:
            pass
        
        elif len(i) == 2:

            if mergetest(i[0],i[1]) == True:
                mergenodes(i[0],i[1])

        else:
            #Create copy of list of node names on the level
            nlist = list(i)
            
            #n is the number of nodes on that level
            n = len(nlist)

            #k list of 0s to mark which nodes have been deleted
            k = np.zeros(n)            
            
            #Create empty numpy matrix of similarity scores
            simmatrix = np.zeros((n,n))

            #Calculate values in half of matrix
            for a, b in itertools.combinations(enumerate(i), 2):
                 simmatrix[a[0],b[0]] = mergetest(a[1],b[1])
             
            #While some values in matrix are not equal to zero
            while np.count_nonzero(simmatrix) != 0:

                #Mask 0s to find minimum value which is not 0
                ma = np.ma.masked_equal(simmatrix, 0.0, copy=False)

                #Return 2D coordinates for position of minimum score
                mnodes = np.unravel_index(ma.argmin(), (n,n))

                #Merge nodes with lowest score
                mergenodes(nlist[mnodes[0]],nlist[mnodes[1]])

                #deleted nodes need to be marked
                k[mnodes[1]] = 1

                #Set all matrix positions of deleted node to 0
                for i in range(len(k)):
                    if i <= mnodes[1]:
                        simmatrix[i][mnodes[1]] = 0
                    elif i > mnodes[1]:
                        simmatrix[mnodes[1]][i] = 0
                                      
                #Recalculate similarity score for changed node
                for i in range(len(k)):
                    if k[i] == 1:
                        continue
                    else:
                        if i == mnodes[0]:
                            continue
                        elif i < mnodes[0]:
                            simmatrix[i][mnodes[0]] = mergetest(nlist[i],nlist[mnodes[0]])
                        elif i > mnodes[0]:
                            simmatrix[mnodes[0]][i] = mergetest(nlist[mnodes[0]],nlist[i])   

    return G







    #sys.stdout.write("Start\t"+str(time.clock())+"\n")

    #parser = argparse.ArgumentParser()

    #parser.add_argument('-t','--tree',dest='tree_input',help='input file for treebuilder')

    #args.tree_input

if 0 != None:
    #Create empty list of allele frequencies
    allelefreq=[]
    alleles=[]

    #Extract genotypes from data
    f = open('/Users/mp18/Documents/bloch/Data/genotype_2.txt', 'rb')
    imputed = []
    for line in f:
        dic = {}   
        l,r= line.split('\t',1)
        imputed.append(l+"\t")
        for i in r:
            if i in ('|','.','/','\t','\n'):
                continue
            else:
                if i in dic:
                    dic[i] += 1
                else:
                    dic[i] = 1
        alleles.append(dic.keys())
        allelefreq.append(dic.values())
    
    #Move current position to beginning of the file
    f.seek(0,0)
    reader = csv.reader(f, delimiter='\t',skipinitialspace = True)        
    #Create list of tuples which contain genotypes of each individual
    GT = zip(*reader)
    f.close()
    sys.stdout.write("Files read\t"+str(time.clock())+"\n")
    #Remove first tuple of position names
    del GT[0]

    #Haplotype length
    hlength=len(GT[0])-1

print "Files processed\t"+str(time.clock())+"\n"
#sys.stdout.write("Files processed\t"+str(time.clock())+"\n")

print alleles
print allelefreq

G = treebuild(GT)

