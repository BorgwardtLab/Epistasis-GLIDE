*** GLIDE v0.1 ***
Tony Kam-Thong
tony@mpipsykl.mpg.de

(C) Copyright 2011, Tony Kam-Thong [tony@mpipsykl.mpg.de]
 
This program is free software: you can redistribute it and/or modify 
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3.0 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

---------------------------------------------------------------------
SUMMARY

This program is a C and CUDA implementation of tabulating linear regression for an
exhaustive pairwise interaction search on a CUDA enabled GPU. 


For a detailed description of the setting and the search algorithm, please refer to the associated paper:
T. Kam-Thong, C.-A. Azencott, L. Cayton, B. Puetz, A. Altmann, N. Karbalai, P. G. Saemann, B. Schoelkopf, B. Mueller-Myhsok, and K. M. Borgwardt. GLIDE: GPU-Based Linear Regression for Detection of Epistasis. Human Heredity, 73 (4) 2012.

---------------------------------------------------------------------
FILES
    Makefile				Makefile

    GLIDE_main.cu			GLIDE main code
    GLIDE_kernel_op.cu			GLIDE kernel code
    CUDAbook.h				CUDA header file
    CUDAcheck.h				CUDA header file
    GLIDE.h				GLIDE header file

    search_name_onestudy.R		R file for post-processing of GLIDE output

    plink2glide.py			Python file for transforming binary PLINK data into a genotype file readable by GLIDE  

    Test1kind_first1ksnp.txt		Test file with 1000 SNPs, 1000 individuals
    Test1kind_second1ksnp.txt 		Test file with 1000 SNPs, 1000 individuals
    Test1kind_pheno.txt			Test file with 1000 phenotypes
    first1k_snpnames.txt		Names of the SNPs in Test1kind_first1ksnp.txt
    second1k_snpnames.txt		Names of the SNPs in Test1kind_second1ksnp.txt
---------------------------------------------------------------------



USAGE

*********************************************
* Make                                      *
*********************************************
make all
(Note:If running on an older GPU of compute capability 1.3,
Please replace -sm_20 to -sm_13 in the Makefile and replace BlockSize 16 to 10 in test.h)

*********************************************
* Calling GLIDE	            		    *
*********************************************

./GLIDE -f1 genoFile -f2 genoFile2 -fp phenoFile -n NSubj -m NSNPS -m2 NSNPS2 -p NSNPS_GPU -t t_thres -o results -g deviceOrdinal

genoFile      = txt file containing the first genotype file in row major order (i.e. each SNP is a row and each subject is a column)\n");
genoFile2     = txt file containing the second genotype file in row major order (i.e. each SNP is a row and each subject is a column)\n");
phenoFile     = txt file containing the phenotype\n");
NSUB	      = Number of Subjects\n");
NSNPS	      = Number of SNPs\n");
NSNPS_GPU     = Size of the partition SNPs, must be integer multiple of BlockSize\n");
t_thres	      = t-score threshold\n");
results       = output file ; stored in text format\n");
deviceNum     = GPU ID # use for the run\n");

NOTE:
1-Selection of the t-threshold (-t) is based on the constraint imposed by the
available storage space on the host machine.  
2-Selection of the partition size (-p) is based on the constraint imposed by the
available memory on the GPU.  
Example-retaining 1 million pairs resulted in file size of 56 MB and for a 1GB of device memory, a partition size of 2000
SNPs have worked reliably well without incurring segmentation fault.          


***************************************************
* Post Processing: assigning snpnames and pvalues *
***************************************************
R --vanilla --args "Resultsname" "Set1_snpname" "Set2_snpname" "chunksize1" "chunksize2" "nsubjects" "NewResultsname" "BlockSize"  < search_name_onestudy.R


---------------------------------------------------------------------
SAMPLE CALL

make all
./GLIDE -f1 Test1kind_first1ksnp.txt -f2 Test1kind_second1ksnp.txt -fp Test1kind_pheno.txt -n 1000 -m 1000 -m2 1000 -p 1000 -t 4 -o Results_1k.txt -g 0
R --vanilla --args "Results_1k.txt" "first1k_snpnames.txt" "second1k_snpnames.txt" "1000" "1000" "1000" "Results_snpnames.txt" "16" < search_name_onestudy.R

---------------------------------------------------------------------
SAMPLE RESULTS
P1 P2 bidx bidy tidx tidy Tint TSnp1 TSnp2 TSnp1n2 Snp1 Snp2 Pint PSnp1 PSnp2 PSnp1n2
0 0 2 47 11 15 1.07445 -1.70963 -4.05737 4.19323 Set1rs44 Set2rs768 0.282881265084217 0.0876457911241603 5.35104201455988e-05 2.99492307327786e-05
0 0 3 31 11 14 1.56592 -3.40932 -3.44837 4.50097 Set1rs60 Set2rs511 0.117684988548390 0.000677358484171921 0.00058758464263497 7.56251433068282e-06
0 0 5 12 3 2 0.86708 -2.55246 -2.88948 4.01919 Set1rs84 Set2rs195 0.386107019174172 0.0108447652338024 0.00394246680550002 6.28008695821487e-05
...


Where:
P1=Partition 1
P2=Partition 2
bidx=block ID x
bidy=block ID y
tidx=thread ID x
tidy=thread ID y
Tint=Tscore of estimated Intercept 
TSnp1=Tscore of estimated SNP1 coefficient
TSnp2=Tscore of estimated SNP2 coefficient
TSnp1n2=Tscore of estimated Interaction coefficient
Snp1= SNP1 name
SNP2= SNP2 name
Pint=p-value of estimated Intercept 
PSnp1=p-value of estimated SNP1 coefficient
PSnp2=p-value of estimated SNP2 coefficient
PSnp1n2=p-value of estimated Interaction coefficient
