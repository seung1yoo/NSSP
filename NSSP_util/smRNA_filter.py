
def main(args):
    out = open(args.outFn, 'w')
    for line in open(args.inFn):
        items = line.rstrip('\n').split('\t')
        if items[0] in ['smRNA_candidate']:
            title = line
            out.write(title)
            idxDic = dict()
            for idx, item in enumerate(items):
                if item.startswith('ReadCnt:'):
                    idxDic.setdefault(idx, item)
            continue
        depths = [int(items[didx]) for didx in idxDic.keys()]
        total_depth = sum(depths)
        if total_depth >= args.depth:
            out.write('{0}\n'.format('\t'.join(items)))

    out.close()

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-ifn', '--inFn',
            default='../do_mergeExpGene/smRNAs.xls')
    parser.add_argument('-d', '--depth', type = int,
            default=100)
    parser.add_argument('-ofn', '--outFn',
            default='../do_mergeExpGene/smRNAs.d100.xls')
    args = parser.parse_args()
    main(args)
