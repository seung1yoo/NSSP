#!/usr/bin/python
import sys
from Bio import SeqIO

def main(args) :
    out = open(args.output, 'w')
    for record in SeqIO.parse(open(args.conseq), 'fasta'):
        seq_ns = record.description.split()
        seq_n = '{0}_{1}'.format(seq_ns[0], seq_ns[1])
        record.id = seq_n
        record.description = ''
        #
        if 'N' in [''.join(list(set([x for x in record.seq])))]:
            #print record.seq
            continue
        if len(record.seq) < 17:
            continue
        if len(record.seq) > 300:
            continue
        SeqIO.write(record, out, 'fasta')
    out.close()

if __name__ == "__main__" :
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('conseq')
    parser.add_argument('output')
    args = parser.parse_args()
    main(args)

