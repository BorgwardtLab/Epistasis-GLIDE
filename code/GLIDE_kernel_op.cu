/* This file is part of the EpiGPULin (EGL) library.
* (C) Copyright 2011, Tony Kam-Thong [tony@mpipsykl.mpg.de]
*/
#ifndef _GLIDE_KERNEL_H_
#define _GLIDE_KERNEL_H_
#define SQ(i,j) ( (i)*(i) )



#include <stdio.h>
#include <math.h>
#include "GLIDE.h"
//#include "driver.h"

__global__ void GLIDE(float* SNPmat, float* SNPmat2, float* Pheno, float* output, int NSNPS_GPU, int NSUB)
{

//Declarations and abbreviations
unsigned int tidx=threadIdx.x;
unsigned int tidy=threadIdx.y;
unsigned int bidx=blockIdx.x;
unsigned int bidy=blockIdx.y;

unsigned int q=BSx;//assuming same size on the dimension BSx=BSy.
__shared__ float A[p][(4*BSx+2)];//Columns of SNPs and SNPs^2
__shared__ float T[(4*BSx+2)][(4*BSx+2)];//Contains all the products needed to reconstruct individual XTXs



//Initializing T (Perhaps use cudaMemset2D)
int x, y;
for(x = 0; x < (4*BSx+2); x ++) {
    for(y = 0; y < (4*BSx+2); y ++) T[x][y] = 0.f;
}

__syncthreads();

unsigned int intm=NSUB/p;  //number of integer multiples of subjects of size blockx*blocky
unsigned int resm=NSUB%p; //the rest of the non multiple 
//unsigned int counts=0; //keep the number of counts in the first for loop
unsigned int endrowcounter1;
unsigned int endrowcounter2;

if(resm==0)
{endrowcounter1=intm;}
else
{endrowcounter1=intm+1;}

////////////////////////////////////////////////////////////////////////////////
//Constructing the A and T matrices in partitions of p individuals
////////////////////////////////////////////////////////////////////////////////
for (int rowcounter = 0; 
	rowcounter < endrowcounter1; 
	rowcounter += 1)
{
if(rowcounter<intm)
{endrowcounter2=p;}
else
{endrowcounter2=resm;}
 for (int indexcounter=0; indexcounter < endrowcounter2;indexcounter +=1)
 {
if(tidy==0)
{
 A[indexcounter][tidx]=SNPmat[IDX(bidx*q+tidx,rowcounter*p+indexcounter,NSUB)];
 A[indexcounter][tidx+BSx+BSy]=pow(A[indexcounter][tidx],2);//square terms
}
if(tidx==0)
{
 A[indexcounter][tidy+BSx]=SNPmat2[IDX(bidy*q+tidy,rowcounter*p+indexcounter,NSUB)];
 A[indexcounter][tidy+2*BSx+BSy]=pow(A[indexcounter][tidy+BSx],2);//square terms
}
if(tidx==0 & tidy==0)
{
 A[indexcounter][2*BSx+2*BSy]=1;//ones for the intercept
 A[indexcounter][2*BSx+2*BSy+1]=Pheno[rowcounter*p+indexcounter];//Phenotype
}
 }
__syncthreads();

////////////////////////////////////////////////////////////////////////////////
//Doing the Matrix Multiplication and Storing the overall correlation matrix in T
////////////////////////////////////////////////////////////////////////////////
 for (int rowcounter2 = 0; rowcounter2  < endrowcounter2; ++rowcounter2)
	{
//First Row
 T[tidx][tidy] += A[rowcounter2][tidx]*A[rowcounter2][tidy];//
 T[tidx][tidy+BSy] += A[rowcounter2][tidx]*A[rowcounter2][tidy+BSy];//
 T[tidx][tidy+3*BSy] += A[rowcounter2][tidx]*A[rowcounter2][tidy+3*BSy];//
 if(tidy==0)
 {
 T[tidx][4*BSy] += A[rowcounter2][tidx];//
 T[tidx][4*BSy+1] += A[rowcounter2][tidx]*A[rowcounter2][4*BSx+1];//
 T[tidx+BSx][4*BSy] += A[rowcounter2][tidx+BSx];//
 T[tidx+BSx][4*BSy+1] += A[rowcounter2][tidx+BSx]*A[rowcounter2][4*BSx+1];//
 }
//Second Row
 T[tidx+BSx][tidy] += A[rowcounter2][tidx]*A[rowcounter2][tidy+BSy]*A[rowcounter2][4*BSx+1];//SPECIAL CASE:to get AiAj*Pheno//
 T[tidx+BSx][tidy+BSy] += A[rowcounter2][tidx+BSx]*A[rowcounter2][tidy+BSy];//
//Third Row
 T[tidx+2*BSx][tidy+BSy] += A[rowcounter2][tidx+2*BSx]*A[rowcounter2][tidy+BSy];//repeat//
 T[tidx+2*BSx][tidy+3*BSy] += A[rowcounter2][tidx+2*BSx]*A[rowcounter2][tidy+3*BSy];//
//Last two single rows
 if(tidy==0 & tidx==0)
 {
 T[4*BSx][4*BSy+1] += A[rowcounter2][4*BSx+1];//
 }
 __syncthreads();
 }

}
////////////////////////////////////////////////////////////////////////////////
//Constructing  4*4 XTX matrices (BS*BS times)
////////////////////////////////////////////////////////////////////////////////
double XTX[8][5];

XTX[0][0]=NSUB;
XTX[0][1]=T[tidx][4*BSx];
XTX[0][2]=T[tidy+BSx][4*BSx];
XTX[0][3]=T[tidx][tidy+BSx];

XTX[1][0]=XTX[0][1];
XTX[1][1]=T[tidx][tidx];
XTX[1][2]=T[tidx][tidy+BSx];
XTX[1][3]=T[tidx+2*BSx][tidy+BSx];

XTX[2][0]=XTX[0][2];
XTX[2][1]=XTX[1][2];
XTX[2][2]=T[tidy+BSx][tidy+BSx];
XTX[2][3]=T[tidx][tidy+3*BSx];

XTX[3][0]=XTX[0][3];
XTX[3][1]=XTX[1][3];
XTX[3][2]=XTX[2][3];
XTX[3][3]=T[tidx+2*BSx][tidy+3*BSx];

////////////////////////////////////////////////////////////////////////////////
//Multiplication XTY
////////////////////////////////////////////////////////////////////////////////
XTX[0][4]=T[4*BSx][4*BSy+1];
XTX[1][4]=T[tidx][4*BSy+1];
XTX[2][4]=T[tidy+BSx][4*BSy+1];
XTX[3][4]=T[tidx+BSx][tidy];//SPECIAL:AiAjPheno

////////////////////////////////////////////////////////////////////////////////
//Inverse of XTX
////////////////////////////////////////////////////////////////////////////////

double determ=(-XTX[0][0]*XTX[1][1]*XTX[2][2]*XTX[3][3]+XTX[0][0]*XTX[1][1]*pow(XTX[2][3],2)+XTX[0][0]*SQ(XTX[1][2],2)*XTX[3][3]-2*XTX[0][0]*XTX[1][2]*XTX[1][3]*XTX[2][3]+XTX[0][0]*SQ(XTX[1][3],2)*XTX[2][2]+SQ(XTX[0][1],2)*XTX[2][2]*XTX[3][3]-SQ(XTX[0][1],2)*SQ(XTX[2][3],2)-2*XTX[0][1]*XTX[1][2]*XTX[0][2]*XTX[3][3]+2*XTX[0][1]*XTX[1][2]*XTX[0][3]*XTX[2][3]+	2*XTX[0][1]*XTX[1][3]*XTX[0][2]*XTX[2][3]-2*XTX[0][1]*XTX[1][3]*XTX[0][3]*XTX[2][2]+XTX[1][1]*SQ(XTX[0][2],2)*XTX[3][3]-2*XTX[0][2]*XTX[1][1]*XTX[0][3]*XTX[2][3]-SQ(XTX[0][2],2)*SQ(XTX[1][3],2)+2*XTX[0][2]*XTX[1][3]*XTX[0][3]*XTX[1][2]+XTX[1][1]*SQ(XTX[0][3],2)*XTX[2][2]-SQ(XTX[0][3],2)*SQ(XTX[1][2],2));//This calculates and stores the determinant.

XTX[4][0] = (-XTX[1][1]*XTX[2][2]*XTX[3][3]+XTX[1][1]*SQ(XTX[2][3],2)+SQ(XTX[1][2],2)*XTX[3][3]-2*XTX[1][2]*XTX[1][3]*XTX[2][3]+SQ(XTX[1][3],2)*XTX[2][2])/determ;
XTX[4][1] = -(-XTX[0][1]*XTX[2][2]*XTX[3][3]+XTX[0][1]*SQ(XTX[2][3],2)+XTX[1][2]*XTX[0][2]*XTX[3][3]-XTX[1][2]*XTX[0][3]*XTX[2][3]-XTX[1][3]*XTX[0][2]*XTX[2][3]+XTX[1][3]*XTX[0][3]*XTX[2][2])/determ;
XTX[4][2] =  -(XTX[0][1]*XTX[1][2]*XTX[3][3]-XTX[0][1]*XTX[1][3]*XTX[2][3]-XTX[1][1]*XTX[0][2]*XTX[3][3]+XTX[1][1]*XTX[0][3]*XTX[2][3]+XTX[0][2]*SQ(XTX[1][3],2)-XTX[1][3]*XTX[0][3]*XTX[1][2])/determ;
XTX[4][3] = -(-XTX[0][1]*XTX[1][2]*XTX[2][3]+XTX[0][1]*XTX[1][3]*XTX[2][2]+XTX[1][1]*XTX[0][2]*XTX[2][3]-XTX[1][1]*XTX[0][3]*XTX[2][2]-XTX[1][2]*XTX[0][2]*XTX[1][3]+XTX[0][3]*SQ(XTX[1][2],2))/determ;
XTX[5][0] = -(-XTX[0][1]*XTX[2][2]*XTX[3][3]+XTX[0][1]*SQ(XTX[2][3],2)+XTX[1][2]*XTX[0][2]*XTX[3][3]-XTX[1][2]*XTX[0][3]*XTX[2][3]-XTX[1][3]*XTX[0][2]*XTX[2][3]+XTX[1][3]*XTX[0][3]*XTX[2][2])/determ;
XTX[5][1] = (-XTX[0][0]*XTX[2][2]*XTX[3][3]+XTX[0][0]*SQ(XTX[2][3],2)+SQ(XTX[0][2],2)*XTX[3][3]-2*XTX[0][2]*XTX[0][3]*XTX[2][3]+SQ(XTX[0][3],2)*XTX[2][2])/determ;
XTX[5][2] = -(-XTX[0][0]*XTX[1][2]*XTX[3][3]+XTX[0][0]*XTX[1][3]*XTX[2][3]+XTX[0][2]*XTX[0][1]*XTX[3][3]-XTX[0][2]*XTX[0][3]*XTX[1][3]-XTX[0][3]*XTX[0][1]*XTX[2][3]+SQ(XTX[0][3],2)*XTX[1][2])/determ;
XTX[5][3] = -(XTX[0][0]*XTX[1][2]*XTX[2][3]-XTX[0][0]*XTX[1][3]*XTX[2][2]-XTX[0][2]*XTX[0][1]*XTX[2][3]+SQ(XTX[0][2],2)*XTX[1][3]+XTX[0][3]*XTX[0][1]*XTX[2][2]-XTX[0][3]*XTX[0][2]*XTX[1][2])/determ;
XTX[6][0] = -(XTX[0][1]*XTX[1][2]*XTX[3][3]-XTX[0][1]*XTX[1][3]*XTX[2][3]-XTX[1][1]*XTX[0][2]*XTX[3][3]+XTX[1][1]*XTX[0][3]*XTX[2][3]+XTX[0][2]*SQ(XTX[1][3],2)-XTX[1][3]*XTX[0][3]*XTX[1][2])/determ;
XTX[6][1] = -(-XTX[0][0]*XTX[1][2]*XTX[3][3]+XTX[0][0]*XTX[1][3]*XTX[2][3]+XTX[0][2]*XTX[0][1]*XTX[3][3]-XTX[0][2]*XTX[0][3]*XTX[1][3]-XTX[0][3]*XTX[0][1]*XTX[2][3]+SQ(XTX[0][3],2)*XTX[1][2])/determ;
XTX[6][2] = -(XTX[0][0]*XTX[1][1]*XTX[3][3]-XTX[0][0]*SQ(XTX[1][3],2)-SQ(XTX[0][1],2)*XTX[3][3]+2*XTX[0][1]*XTX[0][3]*XTX[1][3]-SQ(XTX[0][3],2)*XTX[1][1])/determ;
XTX[6][3] = (XTX[0][0]*XTX[1][1]*XTX[2][3]-XTX[0][0]*XTX[1][3]*XTX[1][2]-SQ(XTX[0][1],2)*XTX[2][3]+XTX[0][1]*XTX[0][2]*XTX[1][3]+XTX[0][3]*XTX[0][1]*XTX[1][2]-XTX[0][3]*XTX[0][2]*XTX[1][1])/determ;
XTX[7][0] = -(-XTX[0][1]*XTX[1][2]*XTX[2][3]+XTX[0][1]*XTX[1][3]*XTX[2][2]+XTX[1][1]*XTX[0][2]*XTX[2][3]-XTX[1][1]*XTX[0][3]*XTX[2][2]-XTX[1][2]*XTX[0][2]*XTX[1][3]+XTX[0][3]*SQ(XTX[1][2],2))/determ;
XTX[7][1] = -(XTX[0][0]*XTX[1][2]*XTX[2][3]-XTX[0][0]*XTX[1][3]*XTX[2][2]-XTX[0][2]*XTX[0][1]*XTX[2][3]+SQ(XTX[0][2],2)*XTX[1][3]+XTX[0][3]*XTX[0][1]*XTX[2][2]-XTX[0][3]*XTX[0][2]*XTX[1][2])/determ;
XTX[7][2] = (XTX[0][0]*XTX[1][1]*XTX[2][3]-XTX[0][0]*XTX[1][3]*XTX[1][2]-SQ(XTX[0][1],2)*XTX[2][3]+XTX[0][1]*XTX[0][2]*XTX[1][3]+XTX[0][3]*XTX[0][1]*XTX[1][2]-XTX[0][3]*XTX[0][2]*XTX[1][1])/determ;
XTX[7][3] = -(XTX[0][0]*XTX[1][1]*XTX[2][2]-XTX[0][0]*SQ(XTX[1][2],2)-SQ(XTX[0][1],2)*XTX[2][2]+2*XTX[0][1]*XTX[0][2]*XTX[1][2]-SQ(XTX[0][2],2)*XTX[1][1])/determ;


for (int x = NCOL_OUTPUT_GPU-4; x < NCOL_OUTPUT_GPU;++x)
{
	XTX[x][4]=0.f;

}
////////////////////////////////////////////////////////////////////////////////
//Multiplication invXTXXTY 
////////////////////////////////////////////////////////////////////////////////
for (int colcounter = 0; colcounter < 4; ++colcounter)
{
	XTX[4][4] += XTX[4][colcounter]*XTX[colcounter][4];
	XTX[5][4] += XTX[5][colcounter]*XTX[colcounter][4];
	XTX[6][4] += XTX[6][colcounter]*XTX[colcounter][4];
	XTX[7][4] += XTX[7][colcounter]*XTX[colcounter][4];
}

////////////////////////////////////////////////////////////////////////////////
//Caculating SSE
////////////////////////////////////////////////////////////////////////////////
float SSE=0.f;

for (int rowcounter = 0; 
	rowcounter < endrowcounter1; 
	rowcounter += 1)
{
if(rowcounter<intm)
	{endrowcounter2=p;}
else
	{endrowcounter2=resm;}

 for (int indexcounter=0;
	 indexcounter < endrowcounter2;
	 indexcounter +=1)
 {
	 SSE+=
	  SQ((XTX[4][4]+SNPmat[IDX(bidx*q+tidx,rowcounter*p+indexcounter,NSUB)]*XTX[5][4]+ 
	  SNPmat2[IDX(bidy*q+tidy,rowcounter*p+indexcounter,NSUB)]*XTX[6][4]+
	  SNPmat[IDX(bidx*q+tidx,rowcounter*p+indexcounter,NSUB)]*SNPmat2[IDX(bidy*q+tidy,rowcounter*p+indexcounter,NSUB)]*XTX[7][4]
	  )-Pheno[rowcounter*p+indexcounter],2);
	 __syncthreads();
 }
}
////////////////////////////////////////////////////////////////////////////////
//Storing t-scores back on global memory
////////////////////////////////////////////////////////////////////////////////
unsigned int index = (bidx*((NSNPS_GPU/BSy)*BSx*BSy*NCOL_OUTPUT_GPU) + bidy*(BSx*BSy*NCOL_OUTPUT_GPU))+(tidx*(BSy*NCOL_OUTPUT_GPU)+tidy*NCOL_OUTPUT_GPU);
output[index]=bidx;
output[index +1]=bidy;
output[index +2]=tidx;
output[index +3]=tidy;


output[index +4]=XTX[4][4]/sqrt(XTX[4][0]*SSE/(NSUB-4));
output[index +5]=XTX[5][4]/sqrt(XTX[5][1]*SSE/(NSUB-4));
output[index +6]=XTX[6][4]/sqrt(XTX[6][2]*SSE/(NSUB-4));
output[index +7]=XTX[7][4]/sqrt(XTX[7][3]*SSE/(NSUB-4));


}
#endif
