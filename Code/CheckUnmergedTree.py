import networkx as nx


def CheckUnmergedGraph(G):
    
    #Check that the graph is a directed acyclic graph (DAG)
    if is_directed_acyclic_graph(G) == False:
        print 'Graph is not a directed acyclic graph'
        return False
    else: return True

    #Check that graph is aperiodic
    if is_aperiodic(G) == False:
        print 'Graph is periodic'
        return False
    else: return True

    
    
