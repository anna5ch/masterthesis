#!/usr/bin/env python3
# Christopher Vollmers
# Roger Volden

import sys
import argparse

parser = argparse.ArgumentParser()


parser.add_argument('--infile', '-i', type=str, action='store')
parser.add_argument('--outfile', '-o', type=str, action='store')
parser.add_argument('--primary', '-p', action='store_true')

minimum_intron_size=10


args = parser.parse_args()
psl_file = args.infile
clean_psl_file = args.outfile
primary=args.primary


def parse_contigs(psl_file, clean_psl_file,primary):
    '''
    Adjusts some of the psl values
    '''
    used=set()
    out = open(clean_psl_file, 'w')
    for line in open(psl_file):
        a = line.strip().split('\t')
        start = int(a[15])
        blocksizes = a[18].split(',')[:-1]
        blockstarts = a[20].split(',')[:-1]

        name=a[9]
        if primary and name in used:
            continue

        blockstarts_clean = []
        readstarts_clean = []
        blocksizes_clean = []
        size_gap = []

        for x in range(len(blocksizes)):
            blockstart = int(blockstarts[x])
            blocksize = int(blocksizes[x])
            blockend = blockstart + blocksize
            size_gap.append(blocksize)

            try:
                next_blockstart = int(blockstarts[x + 1])
                gap = next_blockstart - blockend
                size_gap.append(gap)
            except IndexError:
                pass

        new_size_gap = []
        block = 0
        for index in range(len(size_gap)):
            if index % 2 == 0:
                block += size_gap[index]
            if index % 2 == 1:
                if size_gap[index] < minimum_intron_size:
                    block += size_gap[index]
                else:
                    new_size_gap.append(block)
                    new_size_gap.append(size_gap[index])
                    block = 0

        new_size_gap.append(block)

        current_position = start
        current_readposition = int(a[11])
        for index in range(0, len(new_size_gap), 1):
            inter = new_size_gap[index]
            if index % 2 == 0:
                blockstarts_clean.append(str(current_position))
                blocksizes_clean.append(str(inter))
                readstarts_clean.append(str(current_readposition))

                current_position += inter
                current_readposition += inter

            if index % 2 == 1:
                current_position += inter

        a[17] = str(len(blockstarts_clean))

        blockstarts_clean.append('')
        blocksizes_clean.append('')
        readstarts_clean.append('')

        a[18] = (',').join(blocksizes_clean)
        a[19] = (',').join(readstarts_clean)
        a[20] = (',').join(blockstarts_clean)



        new_line = ('\t').join(a)
        out.write(new_line + '\n')
        used.add(name)


def main():
    parse_contigs(psl_file, clean_psl_file,primary)


main()
