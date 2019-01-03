#!/usr/bin/python

def mkGenDic() :
    fr = open(sys.argv[1], "r")

    dic = dict()

    for line in fr.xreadlines() :
        splitted = line.strip().split("\t")
        genType  = splitted[2]

        if genType == "gene" :
            start   = int(splitted[3])
            end     = int(splitted[4])
            genName = splitted[-1].split(";")[0].split()[1].strip('"')

            dic.setdefault(genName, [start, end])

        else :
            continue

    fr.close()
    return dic

def mksRNADic() :
    fr = open(sys.argv[2], "r")

    dic = dict()

    for line in fr.xreadlines() :
        if line.startswith("candidate") :
            continue

        else :
            splitted = line.strip().split("\t")
            sRNAId   = splitted[0]
            start    = int(splitted[1])
            end      = int(splitted[2])
            length   = int(splitted[3])

            dic.setdefault(sRNAId, [start, end, length])

    fr.close()
    return dic

def main() :
    sRNADic = mksRNADic()
    genDic  = mkGenDic()

    dic = dict()

    for sRNA in sRNADic :

        tmp_dic_5prime = dict()
        tmp_dic_3prime = dict()
        tmp_lst        = list()
        check          = None

        for gen in genDic :
            if genDic[gen][0] > sRNADic[sRNA][0] and genDic[gen][1] < sRNADic[sRNA][1] :
                if tmp_lst :
                    dic[sRNA].append(gen)

                else :
                    dic.setdefault(sRNA, ["complete_external", "0", gen])
                    tmp_lst.append(gen)
                    check = "check"
                
                continue

            elif genDic[gen][0] < sRNADic[sRNA][0] and genDic[gen][1] > sRNADic[sRNA][1] :
                dic.setdefault(sRNA, ["complete_internal", "0", gen])
                check = "check"
                continue

            elif genDic[gen][0] > sRNADic[sRNA][0] and genDic[gen][0] < sRNADic[sRNA][1] :
                if tmp_lst :
                    dic[sRNA].append(gen)
            
                else :
                    dic.setdefault(sRNA, ["5prime_partial_internal", "0", gen])
                    check = "check"
                
                continue

            elif genDic[gen][1] > sRNADic[sRNA][0] and genDic[gen][1] < sRNADic[sRNA][1] :
                if tmp_lst :
                    dic[sRNA].append(gen)
            
                else :
                    dic.setdefault(sRNA, ["3prime_partial_internal", "0", gen])
                    check = "check"
                
                continue

            elif genDic[gen][1] < sRNADic[sRNA][0] :
                gap = sRNADic[sRNA][0] - genDic[gen][1]
                tmp_dic_5prime.setdefault(gap, [sRNA, gen])
            
            elif genDic[gen][0] > sRNADic[sRNA][1] :
                gap = genDic[gen][0] - sRNADic[sRNA][1]
                tmp_dic_3prime.setdefault(gap, [sRNA, gen])

        if check == "check" :
            continue

        else :
            if tmp_dic_5prime :
                if tmp_dic_3prime :
                    min_gap_5prime = min(tmp_dic_5prime)
                    sel_gen_5prime = tmp_dic_5prime[min_gap_5prime][1]

                    min_gap_3prime = min(tmp_dic_3prime)
                    sel_gen_3prime = tmp_dic_3prime[min_gap_3prime][1]

                    dic.setdefault(sRNA, ["intergenic", min_gap_5prime, min_gap_3prime, sel_gen_5prime, sel_gen_3prime])

                if not tmp_dic_3prime :
                    min_gap_5prime = min(tmp_dic_5prime)
                    sel_gen_5prime = tmp_dic_5prime[min_gap_5prime][1]

                    dic.setdefault(sRNA, ["intergenic", min_gap_5prime, "0", sel_gen_5prime, "-"])

    fr = open(sys.argv[2], "r")
    fw = open(sys.argv[3], "w")
    fw.write("candidate sRNA ID\tType\tGap\tGene Name\n")

    for line in fr.xreadlines() :
        if line.startswith("candidate") :
            continue

        else :
            splitted = line.strip().split("\t")
            sRNA     = splitted[0]

            if sRNA in dic :
                if dic[sRNA][0] == "intergenic" :
                    fw.write("%s\t%s\t%s,%s\t5prime:%s,3prime:%s\n" % (sRNA, dic[sRNA][0], dic[sRNA][1], dic[sRNA][2], dic[sRNA][3], dic[sRNA][4]))

                else :
                    fw.write("%s\t%s\t%s\t%s\n" % (sRNA, dic[sRNA][0], dic[sRNA][1], ",".join(dic[sRNA][2:])))
            
            else :
                pass
                #print sRNA


    fr.close()
    fw.close()

if __name__ == "__main__" :
    import sys, pdb

    if len(sys.argv) != 4 :
        print "Usage : %s <GenePrediction_mod.gtf> <addReadCnt_result.xls> <output>" % sys.argv[0]
        sys.exit()

    else :
        main()
