import argparse
import numpy as np
import os

parser = argparse.ArgumentParser()
parser.add_argument('file', type=argparse.FileType('r'))
parser.add_argument('outdir', type=str)
args = parser.parse_args()

# keeping this here in case the docker turns out to work different
bigfile = np.asarray(args.file.readlines())

#time to find the header
headind = 0
while bigfile[headind].find('MOTIF')==-1:
    headind += 1
# so ind is where the header ends and the MOTIFs start
ind1 = headind
counter = 1
for i in range((ind1+1),len(bigfile)):
    if bigfile[i].find('MOTIF')>-1:
        ind2 = i
        inds_to_write = np.append(np.arange(headind),np.arange(ind1,ind2))
        bigfile[ind1] = bigfile[ind1].upper()
        lines_to_write = bigfile[inds_to_write]
        motname = bigfile[ind1].split()[1]
        # Use with statement to ensure the file is properly closed
        with open(os.path.normcase(args.outdir + motname + '.txt'), 'w') as writer:
            writer.writelines(lines_to_write)
        # the start of the new motif is here
        ind1 = ind2
        counter = counter+1

# Handling the last motif
ind2 = i+1
inds_to_write = np.append(np.arange(headind),np.arange(ind1,ind2))
bigfile[ind1] = bigfile[ind1].upper()
lines_to_write = bigfile[inds_to_write]
motname = bigfile[ind1].split()[1]

# Use with statement for the last motif as well
with open(os.path.normcase(args.outdir + motname + '.txt'), 'w') as writer:
    writer.writelines(lines_to_write)