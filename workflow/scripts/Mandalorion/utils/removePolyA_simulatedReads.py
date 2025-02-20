import sys
import mappy

def removePolyA(seq):
    reverse = seq[::-1]
    Astate, Astretch, Vstretch, trimPos , Astart = False, 0, 0, 0, 0
    for pos in range(0, len(reverse), 1):
        base = reverse[pos]
        if not Astate:
            if base == 'A':
                Astretch += 1
                if Astretch == 6:
                    Astate = True
                    lastA = pos
                    Astart = pos
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
    return seqTrim, Astate, Astart,trimPos


def main():
    trim=False
    if len(sys.argv)>3:
        trim=sys.argv[3].split(',')
        fivePrime=int(trim[0])
        threePrime=-int(trim[1])
    out=open(sys.argv[2],'w')

    avgQ=15
    coverage=50
    rawLength=100000

    for name,seq,qual in mappy.fastx_read(sys.argv[1]):
        rootName=name.replace('_','-')
        if trim:
            seq=seq[fivePrime:threePrime]
        length=len(seq)
#        newName=rootName+'_'+str(avgQ)+'_'+str(rawLength)+'_'+str(coverage)+'_'+str(length)+'_'+str(length)
        newName=name
        use=False
        seqTrimF,AstateF, AstartF, trimPosF = removePolyA(seq)
        seqTrimR,AstateR, AstartR, trimPosR = removePolyA(mappy.revcomp(seq))
        if AstateF and not AstateR:
            seqTrim = seqTrimF
            use=True
        elif AstateF and AstateR:
            if AstartF < AstartR:
                seqTrim = seqTrimF
                use=True
            elif AstartR < AstartF:
                seqTrim = seqTrimR
                use = True
            elif AstartR == AstartF:
                if trimPosF>trimPosR:
                    seqTrim=seqTrimF
                    use=True
                elif trimPosR>trimPosF:
                    seqTrim=seqTrimR
                    use=True
        elif AstateR and not AstateF:
            use=True
            seqTrim = seqTrimR
        if use:
            out.write('>'+newName+'\n'+seqTrim+'\n')
#        else:
#            print('>'+newName, AstateF,AstateR,trimPosF,trimPosR)
#            print(seq)
    out.close()

main()
