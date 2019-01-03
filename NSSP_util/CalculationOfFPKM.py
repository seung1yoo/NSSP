import sys

if len(sys.argv) != 4 :
	print "Usage : python %s <.bowtie.sort.idxstats> <.FPKM.value> <read.count>" % sys.argv[0]
	sys.exit()

f = open(sys.argv[1], 'r')

LENGTH = {}
READCNT = {}
TOTALMAPPED = 0

for lines in f.xreadlines() :
	words = lines.rstrip().split("\t")
	if words[0] == '*' :
		continue
	LENGTH[words[0]] = int(words[1])
	READCNT[words[0]] = int(words[2])
	TOTALMAPPED += int(words[2])

f.close()

w1 = open(sys.argv[2], 'w')
w2 = open(sys.argv[3], 'w')
w1.write("#UnigeneID\tFPKM\n")
w2.write("#UnigeneID\tReadCnt\n")

for SEQID in sorted(LENGTH.keys()) :
	FPKMVALUE = (1000000000.0 * (READCNT[SEQID]/2) / ((TOTALMAPPED/2) * LENGTH[SEQID]))
	w1.write("%s\t%.2f\n" % (SEQID, FPKMVALUE))
	w2.write("%s\t%d\n" % (SEQID, READCNT[SEQID]))

w1.close()
w2.close()
