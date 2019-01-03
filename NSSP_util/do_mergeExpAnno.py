
def exp_parser(merge_dic, expFn):
    for line in open(expFn):
        items = line.rstrip('\n').split('\t')
        if items[0] in ['candidate sRNA ID']:
            idxDic = dict()
            samples = list()
            for idx, item in enumerate(items):
                idxDic.setdefault(item, idx)
                if item.startswith('ReadCnt:'):
                    samples.append(item.split(':')[-1])
            continue
        smrna = items[idxDic['candidate sRNA ID']]
        contig = '_'.join(smrna.split('_')[:-1])
        start = items[idxDic['Start']]
        end = items[idxDic['End']]
        length = items[idxDic['Length']]
        conseq = items[idxDic['consensus sequence']]
        exps = [items[idxDic['ReadCnt:{0}'.format(sample)]] for sample in samples]
        #
        merge_dic.setdefault(smrna, {}).setdefault('contig', contig)
        merge_dic.setdefault(smrna, {}).setdefault('start', start)
        merge_dic.setdefault(smrna, {}).setdefault('end', end)
        merge_dic.setdefault(smrna, {}).setdefault('length', length)
        merge_dic.setdefault(smrna, {}).setdefault('conseq', conseq)
        merge_dic.setdefault(smrna, {}).setdefault('samp', samples)
        merge_dic.setdefault(smrna, {}).setdefault('exp', exps)
    return merge_dic, samples

def gene_parser(merge_dic, geneFn):
    for line in open(geneFn):
        items = line.rstrip('\n').split('\t')
        if items[0] in ['candidate sRNA ID']:
            idxDic = dict()
            for idx, item in enumerate(items):
                idxDic.setdefault(item, idx)
            continue
        smrna = items[idxDic['candidate sRNA ID']]
        s_type = items[idxDic['Type']]
        gap = items[idxDic['Gap']]
        gene = items[idxDic['Gene Name']]
        if smrna not in merge_dic:
            print '{0} is not there'.format(smrna)
            import sys
            sys.exit()
        merge_dic[smrna].setdefault('type', s_type)
        merge_dic[smrna].setdefault('gap', gap)
        merge_dic[smrna].setdefault('gene', gene)
    return merge_dic

def makeTable(merge_dic, samples, outdir):
    cols_ = ['smRNA_candidate','contig','start','end','length']
    _cols_ = ['ReadCnt:{0}'.format(sample) for sample in samples]
    _cols = ['type','gap','gene','conseq']
    fn = '{0}/smRNAs.xls'.format(outdir)
    out = open(fn, 'w')
    out.write('{0}\t'.format('\t'.join(cols_)))
    out.write('{0}\t'.format('\t'.join(_cols_)))
    out.write('{0}\n'.format('\t'.join(_cols)))
    for smrna, infoDic in sorted(merge_dic.iteritems()):
        items = [smrna]
        for col in cols_[1:]:
            items.append(infoDic[col])
        items.extend(infoDic['exp'])
        for col in _cols:
            if col in infoDic:
                items.append(infoDic[col])
            else:
                items.append('-')
        out.write('{0}\n'.format('\t'.join(items)))
    out.close()

def main(args):
    merge_dic = dict()
    merge_dic, samples = exp_parser(merge_dic, args.expFn)
    merge_dic = gene_parser(merge_dic, args.geneFn)
    makeTable(merge_dic, samples, args.outdir)

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-efn', '--expFn')
    parser.add_argument('-gfn', '--geneFn')
    parser.add_argument('-od', '--outdir')
    args = parser.parse_args()
    main(args)
