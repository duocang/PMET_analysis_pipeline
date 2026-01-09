import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('promfile',type=str)
parser.add_argument('gff3file',type=str)
parser.add_argument('universefile',type=str)

args = parser.parse_args()

promfile = args.promfile
gff3file = args.gff3file
universefile = args.universefile


# soak up current promoters.bed file
prom = pd.read_csv(promfile,header=None,index_col=None,sep='\t').values
dtype_dict = {0: str, 1: str, 2: str, 3: int, 4: int, 5: str, 6: str, 7: str, 8: str}
annot = pd.read_csv(gff3file, header=None, index_col=None, sep='\t', comment='#', dtype=dtype_dict).values

univ = pd.read_csv(universefile,header=None,index_col=None).values

# the order is the same in the universe.txt as it is in annot.gff3
# after all, annot.gff3 was used to craft universe.txt
# however, promoters.bed might have lost some along the way, so err on the side of caution
geneinds = []
for i in np.arange(annot.shape[0]):
	if 'gene' in annot[i,:]:
		geneinds.append(i)

#and need an "end of file" geneinds as well
geneinds.append(annot.shape[0]+1)

# loop over the genes
for i in np.arange(prom.shape[0]):
	# find which of the geneinds is ours
	j = np.where(univ==prom[i,3])[0][0]
	# dig out the subannot and filter it to cds

	# there are some cases: some part of a genes (several rows) are in other gene's between
	# 1	araport11	gene	589584	590406	.	+	.	ID=gene:AT1G02705;biotype=protei
	# 1	araport11	mRNA	589584	590406	.	+	.	ID=transcript:AT1G02705.1;Parent
	# 1	araport11	exon	589584	590003	.	+	.	Parent=transcript:AT1G02705.1;
	# 1	araport11	CDS	    589584	590003	.	+	0	ID=CDS:AT1G02705.1;Parent=transcript

	# 1	araport11	gene	589706	589996	.	-	.	ID=gene:AT1G02710;biotype=protei
	# 1	araport11	mRNA	589706	589996	.	-	.	ID=transcript:AT1G02710.1;Parent
	# 1	araport11	exon	589706	589996	.	-	.	Parent=transcript:AT1G02710.1;Na
	# 1	araport11	CDS	    589706	589996	.	-	0	ID=CDS:AT1G02710.1;Parent=transcript

	# 1	araport11	exon	590401	590406	.	+	.	Parent=transcript:AT1G02705.1;Na
	# 1	araport11	CDS	    590401	590406	.	+	0	ID=CDS:AT1G02705.1;Parent=transcript
	subannot = annot[geneinds[j]:geneinds[j+1],:]

	noncds = []
	for j in np.arange(subannot.shape[0]):
		if 'CDS' not in subannot[j,:]:
			noncds.append(j)
	temp = np.delete(subannot,noncds,axis=0)

	if temp.size:
		if prom[i,5]=='+':
			# "left to right" transcript. find earliest CDS start and replace end of promoter
			newpos = np.min(temp[:,3])
			prom[i,2] = newpos
		else:
			# "right to left" transcript. find latest CDS "end", which is in fact the earliest CDS start
			# and replace "start" of promoter which is in fact its end
			newpos = np.max(temp[:,4])
			prom[i,1] = newpos

# export the 5' UTR'd promoter file
np.savetxt(promfile, prom, fmt='%s\t%i\t%i\t%s\t%i\t%s')
