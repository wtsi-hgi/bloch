# Copyright (c) 2012 Genome Research Ltd. Author: Michelle Parker <mp18@sanger.ac.uk>
# This file is part of BLOCH.
# BLOCK is free software: you can redistribute it and/or modify it under the terms
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



