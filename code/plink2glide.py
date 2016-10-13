# plink2glide.py -- Convert data in PLINK format to GLIDE input format

# Copyright 2012 Chloe-Agathe Azencott
# Machine Learning and Computational Biology Research Group
# MPI Tuebingen (Germany)


import os
import subprocess
import sys


def main(args):
    usage = """python %s <plink bed root>
    Convert binary PLINK data into input file for GLIDE.
    Produce <plink bed root>.raw (deleted), .snpNames, .pheno, .glideIn
    GLIDE input format: numSNPs lines x numIndividuals columns
    """ % args[0]

    if len(args) != 2:
        print usage
        sys.exit(-1)

    proot = args[1]

    # recode plink data in 0/1/2 format (number of minor alleles)
    # produces <bfile_root>.raw
    cmd = "%s --bfile %s --recodeA --out %s" % (PLINK, proot, proot)    
    returnCode = subprocess.call(cmd, shell=True)
    
    # get the list of SNP names (with '_<minor allele>' appended)
    # one SNP per line
    f = open("%s.raw" % proot)
    hdr = f.readline()
    f.close()

    f = open("%s.snpNames" % proot, 'w')
    f.write("\n".join(hdr.split()[6:]))
    f.close()    
    
    # extract the phenotype
    f = open("%s.raw" % proot)
    g = open("%s.pheno" % proot, 'w')
    hdr = f.readline()
    for line in f:
        g.write("%s\n" % line.split()[5])
    f.close()
    g.close()

    # extract the genotype
    f = open("%s.raw" % proot)
    hdr = f.readline()
    snps = [[] for i in range(len(hdr.split()[6:]))]
    for line in f:
        for (i, snpVal) in enumerate(line.split()[6:]):
            snps[i].append(snpVal)
    f.close()

    g = open("%s.glideIn" % proot, 'w')
    for snpVals in snps:
        g.write("%s\n" % " ".join(snpVals))        
    g.close()
    

    # remove raw file
    os.remove("%s.raw" % proot)



if __name__ == "__main__":
    main(sys.argv)


