def function(a,b):
    if a < b:
        return False
    else:
        return a - b



#Iterate through list and output pairs


levelmin = 100
mini=0
minj=0
for i in range(2,11):
    for j in range(2,11):
        if i != j:
            if function(i,j) < levelmin:
                print 'newminimum'
                print function(i,j)
                print i,j
                levelmin=function(i,j)
                mini=i
                minj=j
                

#Run function which adds the two numbers store minimum.

print 'levelmin'
print levelmin
print 'mini'
print mini
print 'minj'
print minj


