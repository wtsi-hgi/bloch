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
import math


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
            d['nn'] += 1
            G.add_edge(c,d['nn'],a=m,weight=1)
            G.add_node(d['nn'], level=l+1,weight=1)
            gnodes[l+1].append(d['nn'])
            return d['nn']
        else:
            for u,v,k,edata in G.out_edges(c,data=True, keys=True):                
                if edata['a'] == m:
                    G.edge[u][v][k]['weight'] += 1
                    G.node[v]['weight'] +=1
                    return v

            d['nn'] +=1
            G.add_edge(c,d['nn'],a=m,weight=1)
            G.add_node(d['nn'], level=l+1,weight=1)
            gnodes[l+1].append(d['nn'])
            return d['nn']   
    
    def add_genotype(gt):
        G.node[1]['weight'] += 2
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
        #Components of formula are calculated
        ta = G.node[a]['weight']
        tb = G.node[b]['weight']
        thr = math.sqrt(math.pow(ta, -1) + math.pow(tb, -1))
        maxs = 0        
        
        q=[[a,b]]       

        while len(q) != 0:

            da = dict((edata['a'],G.edge[u][v][k]['weight']) for u,v,k,edata in G.out_edges(q[0][0], data=True, keys= True))
            db = dict((edata['a'],G.edge[u][v][k]['weight']) for u,v,k,edata in G.out_edges(q[0][1], data=True, keys = True))    

            sa = set(da.keys())           
            sb = set(db.keys())         

            #Iterate through elements that are in sa that are not in sb
            for i in sa.difference(sb):                     
                s = math.fabs(float(da[i])/float(ta))                
                if s > thr:
                    return 0.0
                elif s > maxs:
                    maxs = s

            #Iterate through elements that are in sb but not sa
            for i in sb.difference(sa):        
                s = math.fabs(float(db[i])/float(tb))                
                if s > thr:
                    return 0.0
                elif s > maxs:
                    maxs = s

            #Iterates through alleles which are a member of both node's outgoing edges
            for i in sa.intersection(sb):
                s = math.fabs(float(da[i])/ta - (float(db[i])/tb))
                
                if s > thr:
                    return 0.0             
                else:
                    if s > maxs:
                        maxs = s                        
                    a = [v for u,v,edata in G.out_edges(q[0][0],data=True) if edata['a']==i]                   
                    b = [v for u,v,edata in G.out_edges(q[0][1],data=True) if edata['a']==i]
                    
                    #Adds next nodes to list if edges has passed edge test            
                    if G.node[a[0]]['level'] != (hlength+1):                        
                        q.append([a[0],b[0]])

            q.pop(0)
        return maxs

            

    def mergenodes(a,b):
        print mergetest(a,b)
        if a == 18:
            print mergetest(19,14419)
        print "Node a"
        for i in G.in_edges(a, data=True):
            print i
        print "Node b"
        for i in G.in_edges(b, data=True):
            print i

        #Move all incoming edges of b to a
        for i,j,k,edata in G.in_edges(b, data=True, keys=True):
            G.add_edge(i,a,a=edata['a'],weight=edata['weight'])
            G.node[a]['weight'] += G.edge[i][j][k]['weight']
            G.remove_edge(i,j,key=k)

        q=[[a,b]]

        while len(q) != 0:

            #Pairwise iterate through outgoing edges of a and b
            for ba, bb, bk, bedata in G.out_edges(q[0][1], data=True, keys=True):

                found = 0
                for aa, ab, ak, aedata in G.out_edges(q[0][0], data=True, keys=True):

                    if bedata['a'] == aedata['a']:
                        found += 1
                        G.edge[aa][ab][ak]['weight'] += bedata['weight']
                        G.node[ab]['weight'] += bedata['weight']
                        if G.node[bb]['level'] != (hlength +1):  
                            q.append([ab,bb])
                            
                        break               

                if found == 0:
                    G.add_edge(q[0][0],bb,a=bedata['a'],weight=bedata['weight'])
                   

                if found > 1:
                    print "ERROR"

                G.remove_edge(ba,bb,key=bk)
                if G.node[bb]['level'] == (hlength +1):
                    gnodes[G.node[bb]['level']].remove(bb)
                    G.remove_node(bb)

            gnodes[G.node[q[0][1]]['level']].remove(q[0][1])
            G.remove_node(q[0][1])
            q.pop(0)
      
    
    gnodes = []
    #Create list of nodes on each level.
    for i in range(hlength+2):
        gnodes.append([])    

    #Create networkx MultiGraph
    G=nx.MultiDiGraph()
    
    #Add start node 1.
    G.add_node(1,level=0,weight=0)    
    gnodes[0].append(1)
    d = {'nn': 1}    
    
    for i in GT:
        add_genotype(i)
    print G.order()
    
        
    for i in list(gnodes[:-1]):
        if len(i) == 1:
            pass
        
        elif len(i) == 2:
            if mergetest(i[0],i[1]) != 0.0:
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
                


    d['nn'] += 1
    G.add_node(d['nn'], weight = 0, level=hlength+1)

    for i in list(gnodes[-1]):
        for a,b,k,edata in G.in_edges(i, data=True, keys=True):
            G.add_edge(a,d['nn'],a=edata['a'],weight=edata['weight'])
            G.node[d['nn']]['weight'] += edata['weight']
            G.remove_edge(a,b,key=k)
            gnodes[G.node[i]['level']].remove(i)
            G.remove_node(i)
    gnodes[-1] = [d['nn']]                         
            
    return G

sys.stdout.write("Start\t"+str(time.clock())+"\n")

parser = argparse.ArgumentParser()
parser.add_argument('-t','--tree',dest='tree_input',help='input file for treebuilder')
args = parser.parse_args()   

if args.tree_input != None:
    #Create empty list of allele frequencies
    allelefreq=[]
    alleles=[]

    #Extract genotypes from data
    f = open(args.tree_input, 'rb')
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

sys.stdout.write("Files processed\t"+str(time.clock())+"\n")
G = treebuild(GT)
print G.order()
sys.stdout.write("Tree Built\t"+str(time.clock())+"\n")

