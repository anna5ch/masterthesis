import sys
import numpy as np
import os


inFolder=sys.argv[1]
outFolder=sys.argv[2]
outFolder2=sys.argv[3]
sampleIDs=sys.argv[4].split(',')

pslFile=inFolder+'/Isoforms.filtered.clean.psl'
gtfFile=inFolder+'/Isoforms.filtered.clean.gtf'
quantFile=inFolder+'/Isoforms.filtered.clean.quant'
read2isoforms=inFolder+'/tmp/reads2isoforms.txt'

os.system('scp '+gtfFile+' '+outFolder+'/models.gtf')
os.system('scp '+gtfFile+' '+outFolder2+'/models.gtf')

modelSet=set()
for line in open(pslFile):
    a=line.strip().split('\t')
    name=a[9]
    modelSet.add(name)

readsOut=open(outFolder+'/read_model_map.tsv','w')
readsOut.write('read_id\ttranscript_id\n')
for line in open(read2isoforms):
    a=line.strip().split('\t')
    isoform=a[1]
    if isoform in modelSet:
        readsOut.write(line)
readsOut.close()


inFile=open(quantFile,'r')
first=inFile.readline()

quantOut=open(outFolder2+'/expression.tsv','w')
quantOut.write('ID\t')
for sample in sampleIDs:
    quantOut.write(sample+'\t')
quantOut.write('\n')

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

inFile=open(quantFile,'r')
first=inFile.readline()

while True:
    line=inFile.readline()
    if not line:
        break
    a=line.strip().split('\t')
    isoform=a[0]
    values=a[1:]
    converted_values=[]
    quantOut.write(isoform+'\t')
    for index in range(0,len(values),1):
        converted=int(values[index])*conversionDict[index]
        converted_values.append(str(converted))
    quantOut.write('\t'.join(converted_values)+'\n')

quantOut.close()
