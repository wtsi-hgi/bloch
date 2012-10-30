import networkx as nx
import sys

haplotypes = []
frequency = []

with open('/Users/mp18/Documents/Python/Data/Table_1', 'r') as f:
    data = f.read()
f.close()

count = 0
for word in data.split():
    if count==0:
        num = int(word)
        haplotypes.append(num)
        count = 1
        continue
    if count==1:
        num = int(word)
        frequency.append(num)
        count = 0
        continue

G = nx.MultiDiGraph()



