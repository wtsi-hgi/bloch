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
        for i in range(len(word)):
           
           if  word[i] not in ("1", "2"):
                               
               raise ValueError("Haplotypes contain value that is not 1 or 2")
         
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


    
