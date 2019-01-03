import os
import sys

class NSSP:
    def __init__(self, conf_fn, cpu):
        self.ioput_dic = dict()
        self.conf_dic = self.conf_parser(conf_fn)
        self.cpu = cpu

    def conf_parser(self, conf_fn):
        conf_dic = dict()
        for line in open(conf_fn):
            if line.startswith('#'):
                continue
            items = line.rstrip('\n').split()
            samp_id = items[0]
            samp_fn = items[1]
            if samp_id in conf_dic:
                print 'sample id is duplicated'
                import sys; sys.exit()
            if not os.path.isfile(samp_fn):
                print 'sample file is not in {0}'.format(samp_fn)
                import sys; sys.exit()
            conf_dic.setdefault(samp_id, samp_fn)
            #
            self.ioput_dic.setdefault('RAW', {}).setdefault(samp_id, samp_fn)
        return conf_dic

    def do_cutadapt(self, samp_id):
        outdir = './do_cutadapt/{0}'.format(samp_id)
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        #
        self.ioput_dic.setdefault('cutadapt', {}).setdefault(samp_id,
                '{0}/{1}.clean.fastq'.format(outdir,samp_id))
        #
        alarm('cutadat', [samp_id], '{0}/{1}.clean.fastq'.format(outdir,samp_id))
        #
        if os.path.isfile(self.ioput_dic['cutadapt'][samp_id]):
            pass
        else:
            cmd = 'cutadapt -O 3 -q 20 -m 17 -b TGGAATTCTCGGGTGCCAAGG' \
                  ' --quality-base 33 -f fastq' \
                  ' --untrimmed-output {0}/{1}.fragments.fastq' \
                  ' --too-short-output {0}/{1}.short.fastq' \
                  ' --output {0}/{1}.clean.pre.fastq' \
                  ' {2} > {0}/{1}.clean.pre.fastq.log' \
                  ''.format(outdir, samp_id, self.ioput_dic['RAW'][samp_id])
            print cmd
            os.system(cmd)
            cmd = 'cutadapt' \
                  ' -u 4 -u -4' \
                  ' --output {0}/{1}.clean.fastq' \
                  ' {0}/{1}.clean.pre.fastq > {0}/{1}.clean.fastq.log' \
                  ''.format(outdir, samp_id, self.ioput_dic['RAW'][samp_id])
            print cmd
            os.system(cmd)

    def do_bowtie2(self, samp_id):
        outdir = './do_bowtie2/{0}'.format(samp_id)
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        #
        self.ioput_dic.setdefault('bowtie2', {}).setdefault(samp_id,
                '{0}/{1}.bowtie2.bam'.format(outdir, samp_id))
        #
        alarm('bowtie2', [samp_id], '{0}/{1}.bowtie2.bam'.format(outdir,samp_id))
        #
        if os.path.isfile(self.ioput_dic['bowtie2'][samp_id]):
            pass
        else:
            cmd = 'bowtie2 -p {0} -x ref.fa -U {1} -S {2}/{3}.bowtie2.sam'\
                  ' 2> {2}/{3}.bowtie2.sam.log'\
                  ''.format(self.cpu, self.ioput_dic['cutadapt'][samp_id], outdir, samp_id)
            print cmd
            os.system(cmd)
            cmd = 'samtools view -Sb -F 2048 {0}/{1}.bowtie2.sam > {0}/{1}.bowtie2.bam'\
                  ''.format(outdir, samp_id)
            print cmd
            os.system(cmd)
            cmd = 'samtools sort -o {0}/{1}.bowtie2.sorted.bam {0}/{1}.bowtie2.bam'\
                  ''.format(outdir, samp_id)
            print cmd
            os.system(cmd)
            cmd = 'samtools stats {0}/{1}.bowtie2.sorted.bam > {0}/{1}.bowtie2.sorted.bam.stats'\
                  ''.format(outdir, samp_id)
            print cmd
            os.system(cmd)

    def do_extConSeq(self, samp_ids):
        outdir = './do_extConSeq'
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        #
        self.ioput_dic.setdefault('extConSeq', {}).setdefault('merged',
                '{0}/merged.consensus_mod.fa'.format(outdir))
        #
        alarm('extConSeq', samp_ids, '{0}/merged.consensus_mod.fa'.format(outdir))
        #
        if os.path.isfile(self.ioput_dic['extConSeq']['merged']):
            pass
        else:
            cmd = 'samtools merge -f {0}/merged.bam {1}'\
                  ''.format(outdir, ' '.join([self.ioput_dic['bowtie2'][samp_id] for samp_id in samp_ids]))
            print cmd
            os.system(cmd)
            cmd = 'samtools sort -o {0}/merged.sorted.bam {0}/merged.bam'\
                  ''.format(outdir)
            print cmd
            os.system(cmd)
            cmd = 'python ./NSSP_util/SAMTOOLSConsensusSequence.v.1.0.1.py'\
                  ' ref.fa'\
                  ' {0}/merged.sorted.bam'\
                  ' {0}/merged.consensus.fa'\
                  ''.format(outdir)
            print cmd
            os.system(cmd)
            cmd = 'python ./NSSP_util/modMergedConsensus.py {0}/merged.consensus.fa {0}/merged.consensus_mod.fa'\
                  ''.format(outdir)
            print cmd
            os.system(cmd)
            cmd = 'bowtie2-build {0}/merged.consensus_mod.fa {0}/merged.consensus_mod.fa'\
                  ' 2>&1 > {0}/merged.consensus_mod.fa.bt2.log'\
                  ''.format(outdir)
            print cmd
            os.system(cmd)

    def do_expConSeq(self, samp_id):
        outdir = './do_expConSeq/{0}'.format(samp_id)
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        #
        self.ioput_dic.setdefault('expConSeq', {}).setdefault(samp_id,
                '{0}/{1}.bowtie2.sorted.bam.idxstats'.format(outdir, samp_id))
        #
        alarm('expConSeq', [samp_id], '{0}/{1}.bowtie2.sorted.bam.idxstats'.format(outdir, samp_id))
        #
        if os.path.isfile(self.ioput_dic['expConSeq'][samp_id]):
            pass
        else:
            cmd = 'bowtie2 -p {0} -x {4} -U {1} -S {2}/{3}.bowtie2.sam'\
                  ' 2> {2}/{3}.bowtie2.sam.log'\
                  ''.format(self.cpu, self.ioput_dic['cutadapt'][samp_id], outdir, samp_id, self.ioput_dic['extConSeq']['merged'])
            print cmd
            os.system(cmd)
            cmd = 'samtools view -Sb -F 2048 {0}/{1}.bowtie2.sam > {0}/{1}.bowtie2.bam'\
                  ''.format(outdir, samp_id)
            print cmd
            os.system(cmd)
            cmd = 'samtools sort -o {0}/{1}.bowtie2.sorted.bam {0}/{1}.bowtie2.bam'\
                  ''.format(outdir, samp_id)
            print cmd
            os.system(cmd)
            cmd = 'samtools stats {0}/{1}.bowtie2.sorted.bam > {0}/{1}.bowtie2.sorted.bam.stats'\
                  ''.format(outdir, samp_id)
            print cmd
            os.system(cmd)
            cmd = 'samtools index {0}/{1}.bowtie2.sorted.bam'\
                  ''.format(outdir, samp_id)
            print cmd
            os.system(cmd)
            cmd = 'samtools idxstats {0}/{1}.bowtie2.sorted.bam > {0}/{1}.bowtie2.sorted.bam.idxstats'\
                  ''.format(outdir, samp_id)
            print cmd
            os.system(cmd)
            cmd = 'python ./NSSP_util/CalculationOfFPKM.py {0}/{1}.bowtie2.sorted.bam.idxstats {0}/{1}.bowtie2.sorted.bam.FPKM.xls {0}/{1}.bowtie2.sorted.bam.ReadCnt.xls'\
                  ''.format(outdir, samp_id)
            print cmd
            os.system(cmd)

    def do_multiExp(self, samp_ids):
        outdir = './do_multiExp'
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        #
        self.ioput_dic.setdefault('multiExp', '{0}/multi.Exp'.format(outdir))
        #
        alarm('multiExp', samp_ids, '{0}/summary.Exp'.format(outdir))
        #
        if os.path.isfile(self.ioput_dic['multiExp']):
            pass
        else:
            cmd = 'python ./NSSP_util/addReadCnt.py -cs {0} -o {1}/multi.Exp -ss {2} -is {3}'\
                  ''.format(self.ioput_dic['extConSeq']['merged'],
                            outdir, ' '.join(samp_ids),
                            ' '.join([self.ioput_dic['expConSeq'][samp_id] for samp_id in samp_ids]))
            print cmd
            os.system(cmd)

    def do_annoGene(self):
        outdir = './do_annoGene'
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        #
        self.ioput_dic.setdefault('annoGene', '{0}/extGenName4sRNA_result.xls'.format(outdir))
        #
        alarm('annoGene', [''], '{0}/extGenName4sRNA_result.xls'.format(outdir))
        #
        if os.path.isfile(self.ioput_dic['annoGene']):
            pass
        else:
            cmd = 'python ./NSSP_util/extGenName4sRNA_multi.py --inFile {0} --gtf ref.gtf -od {1} --refseq ref.fa'\
                  ''.format(self.ioput_dic['multiExp'], outdir)
            print cmd
            os.system(cmd)

    def do_mergeExpAnno(self):
        outdir = './do_mergeExpAnno'
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        #
        self.ioput_dic.setdefault('mergeExpAnno', '{0}/smRNAs.xls'.format(outdir))
        #
        alarm('mergeExpAnno', [''], '{0}/smRNAs.xls'.format(outdir))
        #
        if os.path.isfile(self.ioput_dic['mergeExpAnno']):
            pass
        else:
            cmd = 'python ./NSSP_util/do_mergeExpAnno.py -efn {0} -gfn {1} -od {2}'\
                  ''.format(self.ioput_dic['multiExp'], self.ioput_dic['annoGene'], outdir)
            os.system(cmd)

    def do_summary(self):
        outdir = './do_summary'
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        #
        outdir2 = './do_summary/Raw'
        if not os.path.isdir(outdir2):
            os.makedirs(outdir2)
        self.ioput_dic.setdefault('summary_raw', '{0}/Raw.SeqInfo.Report'.format(outdir2))
        alarm('summary_raw', [''], '{0}/Raw.SeqInfo.Report'.format(outdir2))
        if os.path.isfile(self.ioput_dic['summary_raw']):
            pass
        else:
            fastaStat(self.ioput_dic['RAW'], '{0}/Raw'.format(outdir2))
            linker(self.ioput_dic['RAW'], outdir2)
        #
        outdir2 = './do_summary/Clean'
        if not os.path.isdir(outdir2):
            os.makedirs(outdir2)
        self.ioput_dic.setdefault('summary_clean', '{0}/Clean.SeqInfo.Report'.format(outdir2))
        alarm('summary_clean', [''], '{0}/Clean.SeqInfo.Report'.format(outdir2))
        if os.path.isfile(self.ioput_dic['summary_clean']):
            pass
        else:
            fastaStat(self.ioput_dic['cutadapt'], '{0}/Clean'.format(outdir2))
            linker(self.ioput_dic['cutadapt'], outdir2)
        #
        outdir2 = './do_summary/Mapping'
        if not os.path.isdir(outdir2):
            os.makedirs(outdir2)
        self.ioput_dic.setdefault('summary_mapping', '{0}/Mapping.MapInfo.Report'.format(outdir2))
        alarm('summary_mapping', [''], '{0}/Mapping.MapInfo.Report'.format(outdir2))
        if os.path.isfile(self.ioput_dic['summary_mapping']):
            pass
        else:
            mapStat(self.ioput_dic['bowtie2'], '{0}/Mapping'.format(outdir2))
            linker(self.ioput_dic['bowtie2'], outdir2)
        #
        outdir2 = './do_summary/ExtConSeq'
        if not os.path.isdir(outdir2):
            os.makedirs(outdir2)
        self.ioput_dic.setdefault('summary_extconseq', '{0}/extconseq.SeqInfo.Report'.format(outdir2))
        alarm('summary_extconseq', [''], '{0}/extconseq.SeqInfo.Report'.format(outdir2))
        if os.path.isfile(self.ioput_dic['summary_extconseq']):
            pass
        else:
            fastaStat(self.ioput_dic['extConSeq'], '{0}/extconseq'.format(outdir2))
            linker(self.ioput_dic['extConSeq'], outdir2)
        #
        self.ioput_dic.setdefault('summary_result', '{0}/smRNAs.xls'.format(outdir))
        alarm('summary_result', [''], '{0}/smRNAs.xls'.format(outdir))
        if os.path.isfile(self.ioput_dic['summary_result']):
            pass
        else:
            cmd = 'ln -s {0} {1}'.format(os.path.abspath(self.ioput_dic['mergeExpAnno']), outdir)
            os.system(cmd)

def fastaStat(fastaDic, outprefix):
    files = []
    for sample, fasta in fastaDic.iteritems():
        files.append(fasta)
    cmd = 'python ./NSSP_util/SequenceInfo.fasta.py --outprefix {0} --files {1}'.format(outprefix, ' '.join(files))
    print cmd
    os.system(cmd)

def mapStat(bamDic, outprefix):
    files = []
    for sample, bam in bamDic.iteritems():
        files.append('{0}.sorted.bam.stats'.format('.'.join(bam.split('.bam')[:-1])))
    cmd = 'python ./NSSP_util/bamStatsInfo.py --outprefix {0} --files {1}'.format(outprefix, ' '.join(files))
    os.system(cmd)

def linker(fastaDic, outdir):
    for sample, fasta in fastaDic.iteritems():
        cmd = 'ln -s {0} {1}/{2}'.format(os.path.abspath(fasta), outdir, fasta.split('/')[-1])
        os.system(cmd)


def alarm(pipe, samp_ids, ioFn):
    mssg = '''\

# NSSP : ALARM : START   : {0}
# NSSP : ALARM : SAMPLE  : {1}
# NSSP : ALARM : TAGFILE : {2}
'''.format(pipe, ','.join(samp_ids), ioFn)
    print mssg

def main(args):
    nssp = NSSP(args.conf_fn, args.cpu)
    for pipe in args.pipes:
        if pipe in ['cutadapt']:
            for samp_id in nssp.conf_dic.keys():
                nssp.do_cutadapt(samp_id)
        elif pipe in ['bowtie2']:
            for samp_id in nssp.conf_dic.keys():
                nssp.do_bowtie2(samp_id)
        elif pipe in ['extConSeq']:
            samp_ids = nssp.conf_dic.keys()
            nssp.do_extConSeq(samp_ids)
        elif pipe in ['expConSeq']:
            for samp_id in nssp.conf_dic.keys():
                nssp.do_expConSeq(samp_id)
        elif pipe in ['multiExp']:
            samp_ids = nssp.conf_dic.keys()
            nssp.do_multiExp(samp_ids)
        elif pipe in ['annoGene']:
            nssp.do_annoGene()
        elif pipe in ['mergeExpAnno']:
            nssp.do_mergeExpAnno()
        elif pipe in ['summary']:
            nssp.do_summary()



if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('conf_fn')
    parser.add_argument('--pipes', nargs='+',
            choices=('cutadapt', 'bowtie2','extConSeq','expConSeq','multiExp', 'annoGene', 'mergeExpAnno', 'summary'),
            default=['cutadapt','bowtie2','extConSeq','expConSeq','multiExp', 'annoGene', 'mergeExpAnno', 'summary'])
    parser.add_argument('--cpu', default=10)
    args = parser.parse_args()
    main(args)
