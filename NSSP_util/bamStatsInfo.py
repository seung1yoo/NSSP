
def main(args):
    mapDic = dict()
    for file in args.files:
        for line in open(file):
            items = line.rstrip('\n').split('\t')
            if not items[0] in ['SN']:
                continue
            if items[1] in ['sequences:']:
                total_read = items[2]
            if items[1] in ['reads mapped:']:
                mapped_read = items[2]
            if items[1] in ['reads unmapped:']:
                unmapped_read = items[2]
        mapped_rate = round(float(mapped_read)/int(total_read)*100,2)
        mapDic.setdefault(file,
                [total_read, mapped_read, str(mapped_rate)])
    #
    out = open('{0}.MapInfo.Report'.format(args.outprefix), 'w')
    out.write('Sample\tTotalReads\tMappedReads\tMappedRate\n')
    for file, infos in mapDic.iteritems():
        out.write('{0}\t{1}\n'.format(
            file.split('/')[-1].split('.')[0], '\t'.join(infos)))
    out.close()


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--outprefix')
    parser.add_argument('--files', nargs='+')
    args = parser.parse_args()
    main(args)
