import os
from Bio import SeqIO

def split_addReadCnt_by_contig(inFile, outdir):
    #print(inFile)
    rcfDic = dict()
    rcfhDic = dict()
    for line in open(inFile):
        if line.startswith('candidate'):
            title = line
            continue
        items = line.rstrip('\n').split('\t')
        contig_id = '_'.join(items[0].split('_')[:-1])
        prefix = inFile.split('/')[-1]
        rcfDic.setdefault(contig_id, '{2}/temp/{0}.{1}'.format(prefix, contig_id, outdir))
        rcfhDic.setdefault(contig_id, open('{2}/temp/{0}.{1}'.format(prefix, contig_id, outdir), 'w'))
    #
    #print ('# No.Contig : {0}'.format(len(rcfhDic)))
    #
    for contig_id, fh in rcfhDic.items():
        fh.write(title)
    #
    for line in open(inFile):
        if line.startswith('candidate'):
            title = line
            continue
        items = line.rstrip('\n').split('\t')
        contig_id = '_'.join(items[0].split('_')[:-1])
        rcfhDic[contig_id].write(line)
    #
    for contig_id, fh in rcfhDic.items():
        fh.close()
    #
    return rcfDic

def split_gtf_by_contig(gtf, outdir, gl_dic):
    #print(gtf)
    gfDic = dict()
    gfhDic = dict()
    for line in open(gtf):
        items = line.rstrip('\n').split('\t')
        contig_id = items[0]
        prefix = gtf.split('/')[-1]
        gfDic.setdefault(contig_id, '{2}/temp/{0}.{1}'.format(prefix, contig_id, outdir))
        gfhDic.setdefault(contig_id, open('{2}/temp/{0}.{1}'.format(prefix, contig_id, outdir), 'w'))
    #
    #print ('# No.Contig : {0}'.format(len(gfhDic)))
    #
    for contig_id, fh in gfhDic.items():
        fh.write('{0}\tsiyoo\tgene\t1\t1\t.\t+\t.\tgene_id "GenomeStart"\n'.format(contig_id))
    #
    for line in open(gtf):
        items = line.rstrip('\n').split('\t')
        contig_id = items[0]
        gfhDic[contig_id].write(line)
    #
    for contig_id, fh in gfhDic.items():
        fh.write('{0}\tsiyoo\tgene\t{1}\t{1}\t.\t+\t.\tgene_id "GenomeEnd"\n'.format(contig_id, gl_dic[contig_id]))
        fh.close()
    #
    return gfDic

def exe_extGenName4sRNA(rcfDic, gfDic, outdir):
    results = list()
    contig_ids = list(rcfDic.keys())
    contig_ids.extend(list(gfDic.keys()))
    for contig_id in contig_ids:
        results.append('{1}/extGenName4sRNA_result.{0}'.format(contig_id, outdir))
        if os.path.isfile('{1}/extGenName4sRNA_result.{0}'.format(contig_id, outdir)):
            continue
        #
        if contig_id in rcfDic and contig_id in gfDic:
            rcFile = rcfDic[contig_id]
            gFile = gfDic[contig_id]
            cmd = '/usr/bin/python ./NSSP_util/extGenName4sRNA.py {0} {1} {3}/extGenName4sRNA_result.{2}'.format(
                    gFile, rcFile, contig_id, outdir)
            print(cmd)
            os.system(cmd)
        elif contig_id in rcfDic and not contig_id in gfDic:
            rcFile = rcfDic[contig_id]
            out = open('{1}/extGenName4sRNA_result.{0}'.format(contig_id, outdir), 'w')
            out.write('candidate sRNA ID\tType\tGap\tGene Name\n')
            for line in open(rcFile):
                if line.startswith('candidate'):
                    continue
                items = line.rstrip('\n').split('\t')
                out.write('{0}\tNoGeneInfo\t-\t-\n'.format(items[0]))
            out.close()
        elif not contig_id in rcfDic and contig_id in gfDic:
            out = open('{1}/extGenName4sRNA_result.{0}'.format(contig_id, outdir), 'w')
            out.write('candidate sRNA ID\tType\tGap\tGene Name\n')
            out.close()
        else:
            pass
    return results

def merge_result(results, outdir):
    out = open('{0}/extGenName4sRNA_result.xls'.format(outdir), 'w')
    out.write('candidate sRNA ID\tType\tGap\tGene Name\n')
    for result in results:
        for line in open(result):
            if line.startswith('candidate'):
                continue
            items = line.rstrip('\n').split('\t')
            if items[1] in ['NoGeneInfo']:
                continue
            out.write(line)
    out.close()

def genome_length(refseq):
    gl_dic = dict()
    for record in SeqIO.parse(open(refseq), 'fasta'):
        gl_dic.setdefault(record.id, len(record.seq))
    return gl_dic

def main(args):
    gl_dic = genome_length(args.refseq)
    if not os.path.isdir('{0}/temp'.format(args.outdir)):
        os.system('mkdir {0}/temp'.format(args.outdir))
    rcfDic = split_addReadCnt_by_contig(args.inFile, args.outdir)
    gfDic = split_gtf_by_contig(args.gtf, args.outdir, gl_dic)
    #
    results = exe_extGenName4sRNA(rcfDic, gfDic, args.outdir)
    merge_result(results, args.outdir)

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inFile', default='addReadCnt_result.xls')
    parser.add_argument('-g', '--gtf', default='p3_p13839_Fus_grami_v32.gtf')
    parser.add_argument('-rs', '--refseq', default='ref.fa')
    parser.add_argument('-od', '--outdir')
    args = parser.parse_args()
    main(args)
