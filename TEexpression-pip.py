import sys
import os
import pysam

#==============the input bam file should be sorted by name======================#
# python TEexpression-pip.py <prefix of sorted bam file> <TE annotation file>
# e.g. python TEexpression-pip.py SRRXXXXX B73.structuralTEv2.1.07.2019.filteredTE.subtractexon.plusgenes.chr.sort.gff3
# 2 output files are useful: SRRXXXXX.unique.count.txt (per-element) and SRRXXXXX.TEfamily.counts.txt (per-family)

pre = sys.argv[1] # pre is the prefix of input bam file
ref = sys.argv[2] # modified TE annotation file

bam = pre+'.bam' 
samfile = pysam.AlignmentFile(bam,"rb")

temp = pre+'.converted.bam' # converting original number of reads to 1
outf = pysam.AlignmentFile(temp,"wb",template=samfile)

for line in samfile:
	if 'NH' not in line.tags[-1]:
		line.tags = line.tags
	else:
		realnum = line.tags[-1][-1]
		newtag = ('RH', int(realnum))
		oritag = ('NH', 1)
		line.tags = line.tags[:len(line.tags)-1]
		line.tags = line.tags + [oritag, newtag]
	outf.write(line)

samfile.close()
outf.close()

htsam = pre+'.htseq.sam'
intercounts = 'counts_intermediate_{0}.txt'.format(pre)
os.system('htseq-count -s no -t all -i ID -m union -a 0 -f bam -o {0} {1} {2} > {3}'.format(htsam, temp, ref, intercounts))

filtersam = pre+'.htseq.filtered.sam'
newcount = pre+'.counts.txt'
finalcount = pre+'.TEfamily.counts.txt'

os.system('grep -v "UP\|ambiguous\|__" {0} > {1}'.format(htsam, filtersam))
os.system('perl te_family_mapping_ver6.1.pl {0} {1}'.format(filtersam, pre))
os.system('grep -v "^Gene" {0} | cut -d$"\t" -f1,2,6 > {1}'.format(newcount, finalcount))
os.system('sed -i "1s/.*/TE\tUnique\tTotal/" {0}'.format(finalcount))
os.system('rm -rf *.sam *summary.txt %s'%(newcount))