import sys
import mappy
import argparse


parser = argparse.ArgumentParser()

parser.add_argument('--inFile', '-i', type=str, action='store', help='file to be adapter and poly trimmed')
parser.add_argument('--outFile', '-o', type=str, action='store', help='output file')
parser.add_argument('--trimmedBases', '-t', type=str, action='store', help='in the format [number,number] the number of bases that are going to be trimmed from 5prime and 3prime end (pre-polyA trimming)')

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(0)
args = parser.parse_args()
inFile = args.inFile
outFile = args.outFile
trimmedBases=args.trimmedBases



def removePolyA(seq):
    reverse = seq[::-1]
    Astate, Astretch, Vstretch, trimPos = False, 0, 0, 0
    for pos in range(0, len(reverse), 1):
        base = reverse[pos]
        if not Astate:
            if base == 'A':
                Astretch += 1
                if Astretch == 6:
                    Astate = True
                    lastA = pos
            else:
                Astretch=0
        if Astate:
            if base != 'A':
                Vstretch += 1
                Astretch = 0
            else:
                Astretch += 1
                if Astretch >= 3:
                    Vstretch = 0
                    lastA = pos
            if Vstretch >= 3:
                trimPos = lastA
                break
    reverseTrim = reverse[trimPos:]
    seqTrim = reverseTrim[::-1]
    return seqTrim, Astate


def main():
    trim=False
    if len(sys.argv)>3:
        trim=trimmedBases.split(',')
        fivePrime=int(trim[0])
        threePrime=-int(trim[1])
    out=open(outFile,'w')
    for name,seq,qual in mappy.fastx_read(inFile):
        if trim:
            seq=seq[fivePrime:threePrime]
        seqTrim,Astate=removePolyA(seq)
        out.write('>'+name+'\n'+seqTrim+'\n')
    out.close()

main()
