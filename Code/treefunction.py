import networkx as nx

#These defined before function
h = {'0110': 116, '0000': 21, '0001': 79, '0011': 95, '1011': 152, '1001': 112, '1000': 25}

hlength = 3

#Calculated inside function
hnum = len(h)

#Create list of first node in each level
gl=[1]*(hlength+1)

#Function which takes node and splits node according to haplotype dictionary attached to the node
def nodesplit(G,n,l,level):
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
            m = level[l+1]+1
            #Create edge and label it with first character of key
            G.add_edge(n,m,a=key[0],weight=value)
            #Assign dictionary of suffix of that key to node that has been created.
            G.add_node(m, hap={key[1:]:value})

            #Add node to gl
            for i in range(l+1, hlength+1):
                level[i] += 1


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
nodesplit(G,1,0,gl)

#Merge level 1
merge(1)

print G.nodes(data=True)
for i in G.edges(data=True, keys=True):
    print i
