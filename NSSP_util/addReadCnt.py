#!/usr/bin/python
from Bio import SeqIO
import sys, glob

def mksRNAIdLst(conseq) :
    fr = open(conseq, "r")
    lst = list()
    for line in fr.xreadlines() :
        if line.startswith(">") :
            sRNAId = line.rstrip('\n').lstrip('>')
            lst.append(sRNAId)
    fr.close()
    return lst

def mkReadCntDic(idxstats) :
    Files = idxstats
    dic = dict()
    for File in Files :
        sampId = File.split("/")[-1].split('.')[0]
        fr = open(File, "r")
        for line in fr.xreadlines() :
            splitted = line.strip().split("\t")
            sRNAId   = splitted[0]
            readCnt  = splitted[2]
            dic.setdefault(sampId, {}).setdefault(sRNAId, readCnt)
        fr.close()
    return dic

def _mkPosLenDic(conseq) :
    fr = open(conseq, "r")
    dic = dict()
    for line in fr.xreadlines() :
        if line.startswith(">") :
            sRNAId = line.rstrip('\n').lstrip('>')
            start  = sRNAId.split('_')[-1].split('-')[0]
            end    = sRNAId.split('_')[-1].split('-')[1]
            length = int(end) - int(start) + 1

            dic.setdefault(sRNAId, [start, end, str(length)])
        else :
            seq = line.strip()
            dic[sRNAId].append(seq)
    fr.close()
    return dic

def mkPosLenDic(conseq) :
    fr = open(conseq, "r")
    dic = dict()
    for record in SeqIO.parse(open(conseq), 'fasta'):
        sRNAId = record.id
        start  = sRNAId.split('_')[-1].split('-')[0]
        end    = sRNAId.split('_')[-1].split('-')[1]
        length = int(end) - int(start) + 1
        dic.setdefault(sRNAId, [start, end, str(length)])
        dic[sRNAId].append(str(record.seq))
    return dic

def main(args) :
    sampIdLst  = args.samples
    sRNAIdLst  = mksRNAIdLst(args.extConSeq)
    readCntDic = mkReadCntDic(args.idxstats)
    posLenDic  = mkPosLenDic(args.extConSeq)

    fw = open(args.output, "w")
    fw.write("candidate sRNA ID\tStart\tEnd\tLength\t")
    for sampId in sampIdLst :
        fw.write("ReadCnt:%s\t" % sampId)
    fw.write("consensus sequence\n")
    for sRNAId in sRNAIdLst :
        fw.write("%s\t%s\t" % (sRNAId, "\t".join(posLenDic[sRNAId][:-1])))
        for sampId in sampIdLst :
            fw.write("%s\t" % readCntDic[sampId][sRNAId])
        fw.write("%s\n" % posLenDic[sRNAId][-1])
    fw.close()

if __name__ == "__main__" :
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-cs', '--extConSeq')
    parser.add_argument('-o', '--output')
    parser.add_argument('-ss', '--samples', nargs='+')
    parser.add_argument('-is', '--idxstats', nargs='+')
    args = parser.parse_args()
    main(args)
