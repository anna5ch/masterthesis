import sys
import numpy as np

inFile=open(sys.argv[1],'r')
outFile=open(sys.argv[2],'w')

first=inFile.readline()
outFile.write(first)

quantDict={}
while True:
    line=inFile.readline()
    if not line:
        break
    a=line.strip().split('\t')
    isoform=a[0]
    values=a[1:]
    for index in range(0,len(values),1):
        if index not in quantDict:
            quantDict[index]=[]
        quantDict[index].append(int(values[index]))


inFile.close()
conversionDict={}
for index,values in quantDict.items():
    conversionDict[index]=1000000/sum(values)

print(conversionDict)

inFile=open(sys.argv[1],'r')
first=inFile.readline()

while True:
    line=inFile.readline()
    if not line:
        break
    a=line.strip().split('\t')
    isoform=a[0]
    values=a[1:]
    converted_values=[]
    outFile.write(isoform+'\t')
    for index in range(0,len(values),1):
        converted=int(values[index])*conversionDict[index]
        converted_values.append(str(converted))
    outFile.write('\t'.join(converted_values)+'\n')
