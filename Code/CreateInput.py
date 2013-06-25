from random import randrange

#Create list of haplotypes
haplotype = []

#Create list of haplotype frequencies
frequency = []

#Open data file in read mode
with open('/Users/mp18/Documents/bloch/Data/Table_1', 'r') as f:
    data = f.read()
    
#Close file
f.close()

#Seperate data file line and add each item to corresponding list.
count = 0

for word in data.split():    
    if count==0:
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


l1 = 'pos1'
l2 = 'pos2'
l3 = 'pos3'
l4 = 'pos4'

def allzero(list):
    for i in list:
        if i == 0:
            next
        else:
            return False
    return True

s = 0
while allzero(frequency) == False:
    i = randrange(0,len(frequency))
    if frequency[i] == 0:
        next
    else:
        if s==0:
            l1 += '\t'
            l2 += '\t'
            l3 += '\t'
            l4 += '\t'
            s=1
        elif s==1:
            r = random.randint(0,2)
            if r == 0:
                m = '|'
            else:
                m = '/'
            
            l1 += m
            l2 += m
            l3 += m
            l4 += m            
            s=0

        
        l1 += haplotype[i][0]
        l2 += haplotype[i][1]
        l3 += haplotype[i][2]
        l4 += haplotype[i][3]
        frequency[i] -=1


l1 += '\n'
l2 += '\n'
l3 += '\n'
l4 += '\n'

f = open('/Users/mp18/Documents/bloch/Data/genotype_2.txt', 'w')
f.write(l1)
f.write(l2)
f.write(l3)
f.write(l4)

f.close()
