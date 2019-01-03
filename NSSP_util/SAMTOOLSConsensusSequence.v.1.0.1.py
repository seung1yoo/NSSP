import sys, os, getopt

INBAM = None
OUTFILE = None
REFER = None
DEPTH = 0
GAPSIZE = 100
BEDFILE = ''

#SAMTOOLS = "/BiO/BioTools/samtools/samtools-0.1.19/samtools"
SAMTOOLS = "samtools"


def init():
	global SAMTOOLS, BUFFER_SIZE, INBAM, OUTFILE, REFER, QUALITY, DEPTH, BEDFILE, GAPSIZE
	options, args = getopt.getopt(sys.argv[1:], "s:d:l:")
	for op, p in options:
		print op, p
		if op == "-s": 
			SAMTOOLS = p
		if op == "-d":
			DEPTH = int(p)
		if op == "-l":
			BEDFILE = "-l %s " % p
		if op == "-g":
			GAPSIZE = int(p)

	REFER = args[0]
	INBAM = args[1]
	OUTFILE = args[2]

def countAllele(_ref, _seq):
	alleleDic = {"A":0, "T":0, "G":0, "C":0, "In":0, "Del":0, "N":0, 'InSeq':''}
	JumpCnt = 0
	InSeq = {}

	for base in range(0,len(_seq)):
		if base+JumpCnt >= len(_seq) :
			break
		if _seq[base+JumpCnt] == '.' or _seq[base+JumpCnt] == ',':
			alleleDic[_ref] += 1
		elif _seq[base+JumpCnt] == '^' :
			JumpCnt += 1
		elif _seq[base+JumpCnt] == '$' :
			pass
		elif _seq[base+JumpCnt] == '-' :
			JumpCnt += 1 + int(_seq[base+JumpCnt+1])
		elif _seq[base+JumpCnt] == '*' :
			alleleDic["Del"] += 1
		elif _seq[base+JumpCnt] == '+' :
			try :
				int(_seq[base+JumpCnt+3])
			except :
				try :
					int(_seq[base+JumpCnt+2])
				except :
					Seq = _seq[base+JumpCnt+2:(base+JumpCnt+2)+int(_seq[base+JumpCnt+1])].upper().replace("=",'')
					if not Seq == '' :
						alleleDic["In"] += 1
					JumpCnt += 1 + int(_seq[base+JumpCnt+1])
				else :
					Seq = _seq[base+JumpCnt+3:(base+JumpCnt+3)+int(_seq[base+JumpCnt+1:base+JumpCnt+3])].upper().replace("=",'')
					JumpCnt += 2 + int(_seq[base+JumpCnt+1:base+JumpCnt+3])
					alleleDic["In"] += 1
			else :
				Seq = _seq[base+JumpCnt+4:(base+JumpCnt+4)+int(_seq[base+JumpCnt+1:base+JumpCnt+4])].upper().replace("=",'')
				JumpCnt += 3 + int(_seq[base+JumpCnt+1:base+JumpCnt+4])
				alleleDic["In"] += 1
			if InSeq.has_key(Seq) :
				InSeq[Seq] += 1
			else :
				InSeq[Seq] = 1
		elif _seq[base+JumpCnt] == 'A' or _seq[base+JumpCnt] == 'C' or _seq[base+JumpCnt] == 'G' or _seq[base+JumpCnt] == 'T' :
			alleleDic[_seq[base+JumpCnt]] += 1
		elif _seq[base+JumpCnt] == 'a' :
			alleleDic["A"] += 1
		elif _seq[base+JumpCnt] == 't' :
			alleleDic["T"] += 1
		elif _seq[base+JumpCnt] == 'g' :
			alleleDic["G"] += 1
		elif _seq[base+JumpCnt] == 'c' :
			alleleDic["C"] += 1
		elif _seq[base+JumpCnt] == 'n' or _seq[base+JumpCnt] == 'N' :
			alleleDic["N"] += 1

	for Seq in InSeq.keys() :
		if not Seq == '' :
			if alleleDic["InSeq"] :
				alleleDic["InSeq"] += ",%s:%d" % (Seq, InSeq[Seq])
			else :
				alleleDic["InSeq"] = "%s:%d" % (Seq, InSeq[Seq])
			
	return str(alleleDic["A"]), str(alleleDic["T"]), str(alleleDic["G"]), str(alleleDic["C"]), str(alleleDic["N"]), str(alleleDic["In"]), str(alleleDic["Del"]), alleleDic["InSeq"]

	

if __name__ == "__main__" :
	if len(sys.argv) == 1:
		print "Usage : python %s [-s SAMTOOLS %s] [-d minimum Depth] [-l bed file] [-g maximum gap size] <in.Ref_file> <in.BamFile> <out.fasta>" %(sys.argv[0], SAMTOOLS)
		sys.exit()
	else:
		init()


	result = open(OUTFILE, "w")
	print ("%s mpileup -d 2000000 %s-f %s %s" %(SAMTOOLS, BEDFILE, REFER,  INBAM))
	fi_pileup = os.popen("%s mpileup -d 2000000 %s-f %s %s" %(SAMTOOLS, BEDFILE, REFER, INBAM))


	pre = -1
	position = 0
	prechr = ''
	totaldp = ''
	end = 0
	SeqID = ''
	SeqSeq = ''
	for line in fi_pileup.xreadlines():
		words = line.strip().split("\t")
		if not words[0] == prechr and not prechr == '' :
			result.write("%s-%s legnth=%s cov=%.2f\n%s\n" % (SeqID, str(end), str(len(SeqSeq)), Cov, SeqSeq))
			prechr = words[0]
			SeqID = ''
			SeqSeq = ''
			totaldp = 0
			SeqID = ">%s %s" % (words[0], words[1])
			position = int(words[1])
			pre = -1
			end = 0

		if int(words[1]) < position :
			continue
		position = int(words[1])
		bestAllele = ''
		if len(words) < 5 :
			bestAllele = 'N'
		else :
			seq = words[4]
			quality = words[5]
			a,t,g,c,n,inser,delet,inseq = countAllele(words[2].upper(), seq)
			if (int(a) + int(t) + int(g) + int(c) + int(n) + int(delet) + int(inser)) < DEPTH :
				continue
			total = int(a) + int(t) + int(g) + int(c) + int(n)
			rawtotal = int(a) + int(t) + int(g) + int(c) + int(n) + int(delet) + int(inser)
			bestAllele = ''
			if int(a) > int(t) and int(a) > int(g) and int(a) > int(c) and int(a) > int(n) and int(a) > int(delet) :
				bestAllele = 'A'
			elif int(t) > int(a) and int(t) > int(g) and int(t) > int(c) and int(t) > int(n) and int(t) > int(delet) :
				bestAllele = 'T'
			elif int(g) > int(a) and int(g) > int(t) and int(g) > int(c) and int(g) > int(n) and int(g) > int(delet) :
				bestAllele = 'G'
			elif int(c) > int(a) and int(c) > int(t) and int(c) > int(g) and int(c) > int(n) and int(c) > int(delet) :
				bestAllele = 'C'
			elif int(delet) > int(a) and int(delet) > int(t) and int(delet) > int(g) and int(delet) > int(c) and int(delet) > int(n) :
				bestAllele = ''
			else :
				bestAllele = 'N'
			if inseq :
				Allele = inseq.split(",")
				for preseq in Allele :
					seqs = preseq.split(":")
					if (rawtotal/2.0) < int(seqs[1]) :
						bestAllele += seqs[0]

		if prechr == words[0] and (pre+1) == int(words[1]) :
			SeqSeq += bestAllele
			end = int(words[1])
			pre = int(words[1])
			totaldp += rawtotal
		else :
			if prechr == '' :
				SeqID = ">%s %s" % (words[0], words[1])
				SeqSeq = bestAllele
				end = int(words[1])
				pre = int(words[1])
				prechr = words[0]
				totaldp = rawtotal

			elif not prechr == words[0] :
				Cov = totaldp * 1.0 / len(SeqSeq)
				result.write("%s-%s legnth=%s cov=%.2f\n%s\n" % (SeqID, str(end), str(len(SeqSeq)), Cov, SeqSeq))
				SeqID = ">%s %s" % (words[0], words[1])
				SeqSeq = bestAllele
				end = int(words[1])
				pre = int(words[1])
				prechr = words[0]
				totaldp = rawtotal

			else :
				if end :
					if end > int(words[1]) :
					        gap = int(words[1])-end-1
						SeqSeq += 'N' * gap + bestAllele
						end = int(words[1])
						pre = int(words[1])
					else :
						Cov = totaldp * 1.0 / len(SeqSeq)
						result.write("%s-%s legnth=%s cov=%.2f\n%s\n" % (SeqID, str(end), str(len(SeqSeq)), Cov, SeqSeq))
						SeqID = ">%s %s" % (words[0], words[1])
						SeqSeq = bestAllele
						end = int(words[1])
						pre = int(words[1])
						prechr = words[0]
						totaldp = rawtotal
				else :
					SeqID = ">%s %s" % (words[0], words[1])
					SeqSeq = bestAllele
					end = int(words[1])
					pre = int(words[1])
					prechr = words[0]
					totaldp = rawtotal
			
	Cov = totaldp * 1.0 / len(SeqSeq)
	result.write("%s-%s legnth=%s cov=%.2f\n%s\n" % (SeqID, str(end), str(len(SeqSeq)), Cov, SeqSeq))
	result.close()


