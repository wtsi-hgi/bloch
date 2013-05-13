import networkx as nx

#These defined before function
h = {'0110': 116, '0000': 21, '0001': 79, '0011': 95, '1011': 152, '1001': 112, '1000': 25}

hlength = 3

#Calculated inside function
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


#Create list of first node in each level
gl=[1]*(hlength+1)

#Create networkx MultiGraph
G=nx.MultiDiGraph()

#Add start node 1 and split the first node.
G.add_node(1,hap=h)
nodesplit(1,0)

print G.nodes(data=True)
for i in G.edges(data=True, keys=True):
    print i
