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

import csv
import random

f = open('/Users/mp18/Documents/bloch/Data/genotype_10.txt', 'rb')
g = open('/Users/mp18/Documents/bloch/Data/genotype_10b.txt', 'w')



reader = csv.reader(f, delimiter='\t',skipinitialspace = True)


for line in reader:
    string = ""
    for i in line:
        if "|" not in i:
            string += i
        else:
          
            rand = random.randrange(2)
            if rand == 0:
                a=i[0]
                b=i[2]
            else:
                a=i[2]
                b=i[0]   
                
                
            rand = random.randrange(100)
            if rand == 0:
                a = '.'
                
            rand = random.randrange(100)
            if rand == 0:
                b = '.'
                
            string += ("\t"+a+"/"+b)
    string += "\n"
        
            
    g.write(string)
    
    
g.close
f.close



