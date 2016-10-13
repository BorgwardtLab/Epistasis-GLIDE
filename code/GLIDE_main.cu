/* This file is part of the EpiGPULin (EGL) library.
* (C) Copyright 2011, Tony Kam-Thong [tony@mpipsykl.mpg.de]
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <string.h>
#include <iostream>
#include <assert.h>
#include <cublas.h>
#include <time.h>
#include "CUDAbook.h"
#include "GLIDE.h"
#include "CUDAcheck.h"
//#include "driver.h"


// windows
//#include <windows.h>
//#include <crtdbg.h>


// locals
#include "GLIDE_kernel_op.cu"

void zeroInitfloat(float* data, int size);
void randomInitInt(int* data, int size);
void readfloatData(char *dataFile, unsigned int rows, unsigned int cols, float * data);
void readintData(char *dataFile, unsigned int rows, unsigned int cols, int *data);
void readstringData(char *dataFile, unsigned int rows, unsigned int cols, char *data);
void parseInput(int,char**);//User Settings
int PrintDevices(int deviceCount, int deviceSelected);

int NSNPS;
int NSNPS2;
int NSUB;
int NSNPS_GPU;

float t_thres;
int deviceOrdinal;

char *phenoFile;
char *genoFile;
char *genoFile2;
char *results;

// includes, kernels
////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
//void main(char *genoFile, char *phenoFile)
int main(int argc, char **argv)
{

//Parsing the input parameters
parseInput(argc,argv);


Matrixfloat h_SNPmat;
Matrixfloat h_SNPmat2;
Matrixfloat h_Pheno;
Matrixfloat h_output;
Matrixfloat d_SNPmat;
Matrixfloat d_SNPmat2;
Matrixfloat d_Pheno;
Matrixfloat d_output;

//unsigned int deviceNum=0;
//Check for CUDA machines and choose the one.
cudaError_t err = cudaSuccess;

    int deviceCount = 0;
    err = cudaGetDeviceCount(&deviceCount);
    CheckConditionXR_(err == cudaSuccess, err)

    if (deviceCount == 0) {
        printf("error: no devices supporting CUDA.\n");
        return(-1);
    }

	   
    CheckConditionXR_(deviceOrdinal <= deviceCount, cudaErrorInvalidDevice);

    err = (cudaError_t)PrintDevices(deviceCount, deviceOrdinal);
    CheckConditionXR_(err == cudaSuccess, err)

    err = cudaSetDevice(deviceOrdinal);
    CheckConditionXR_(err == cudaSuccess, err)

printf("*****************\n");
printf("Linear Regression\n");
printf("*****************\n");
fprintf(stderr,"GPU Number %d in use \n",deviceOrdinal);
/*Matrix Memory allocation and calling the CUDA function*/
h_SNPmat.width =NSNPS;
h_SNPmat.height=NSUB;
h_SNPmat2.width =NSNPS2;
h_SNPmat2.height=NSUB;
h_Pheno.height=NSUB;
h_Pheno.width =1;
h_output.height=NSNPS_GPU*NSNPS_GPU*NCOL_OUTPUT_GPU;
h_output.width=1;

//Zero padding
if(NSNPS%NSNPS_GPU>0)
{h_SNPmat.width = (int(NSNPS/NSNPS_GPU)+1)*NSNPS_GPU;}
if(NSNPS2%NSNPS_GPU>0)
{h_SNPmat2.width = (int(NSNPS2/NSNPS_GPU)+1)*NSNPS_GPU;}

//Size declaration
unsigned int h_size_SNPmat  = h_SNPmat.width * h_SNPmat.height;
size_t h_mem_size_SNPmat = sizeof(float) * h_size_SNPmat;

unsigned int h_size_SNPmat2  = h_SNPmat2.width * h_SNPmat2.height;
size_t h_mem_size_SNPmat2 = sizeof(float) * h_size_SNPmat2;

size_t h_mem_size_Pheno= NSUB*sizeof(float);
size_t h_mem_size_output= NSNPS_GPU*NSNPS_GPU*NCOL_OUTPUT_GPU*sizeof(float);

//Size declaration for GPU in Parts
d_SNPmat.width=NSNPS_GPU;d_SNPmat.height=NSUB;
d_Pheno.width=h_Pheno.width;d_Pheno.height=NSUB;
d_output.width=1;d_output.height=NSNPS_GPU*NSNPS_GPU*NCOL_OUTPUT_GPU;

unsigned int d_size_SNPmat  = d_SNPmat.width * d_SNPmat.height;
size_t d_mem_size_SNPmat = sizeof(float) * d_size_SNPmat;
size_t d_mem_size_Pheno= NSUB*sizeof(float);
size_t d_mem_size_output= NSNPS_GPU*NSNPS_GPU*NCOL_OUTPUT_GPU*sizeof(float);


//Allocate host memory 
h_SNPmat.elements = (float*)malloc(h_mem_size_SNPmat); 
h_SNPmat2.elements = (float*)malloc(h_mem_size_SNPmat2); 
//Zero padding
/*memfill(h_SNPmat.elements,h_mem_size_SNPmat,0.f);
memfill(h_SNPmat2.elements,h_mem_size_SNPmat2,0.f);
*/
zeroInitfloat(h_SNPmat.elements,h_SNPmat.width*h_SNPmat.height);
zeroInitfloat(h_SNPmat2.elements,h_SNPmat2.width*h_SNPmat2.height);


h_Pheno.elements = (float*)malloc(h_mem_size_Pheno);

//Read in the matrix
readfloatData(genoFile, NSUB, NSNPS, h_SNPmat.elements);
readfloatData(genoFile2, NSUB, NSNPS2, h_SNPmat2.elements);
readfloatData(phenoFile, NSUB, 1, h_Pheno.elements);
/*SNPNAMES
readstringData(identifier, 1, NSNPS, snpnames1);
readstringData(identifier2, 1, NSNPS, snpnames2);
*/
//randomInitfloat(h_Pheno.elements,h_mem_size_Pheno/sizeof(float));
//randomInitInt(h_SNPmat.elements,h_mem_size_SNPmat/sizeof(int));
fprintf(stderr,"Reading in Matrices into Host Success \n");

//Allocate device memory
cudaMalloc((void**)&d_SNPmat.elements, d_mem_size_SNPmat);
cudaMalloc((void**)&d_SNPmat2.elements, d_mem_size_SNPmat);
cudaMalloc((void**)&d_Pheno.elements, d_mem_size_Pheno);
cudaMalloc((void**)&d_output.elements, d_mem_size_output);

//GPU timing events
cudaEvent_t start, stop;
cudaEventCreate(&start);
cudaEventCreate(&stop);
float elapsedTime;

// Copy host memory to device
cudaMemcpy(d_Pheno.elements, h_Pheno.elements, d_mem_size_Pheno, cudaMemcpyHostToDevice);

//Offset to the host array to shift by NSNPS_GPU at each loop
float* offset1=h_SNPmat.elements; 
float* offset2=h_SNPmat2.elements;

for (int index1count=0;index1count<NSNPS;index1count+=NSNPS_GPU)
{

cudaMemcpy(d_SNPmat.elements, offset1, d_mem_size_SNPmat, cudaMemcpyHostToDevice);

for (int index2count=0;index2count<NSNPS2;index2count+=NSNPS_GPU)
{
cudaMemcpy(d_SNPmat2.elements, offset2, d_mem_size_SNPmat, cudaMemcpyHostToDevice);


//Allocate output vector (done within the loop as it is being cleared)
cudaHostAlloc((void**)&h_output.elements,h_mem_size_output,cudaHostAllocDefault);

// kernel call
dim3 dimBlock(BSx, BSy, BSz );
dim3 dimGrid(d_SNPmat.width/BSx, d_SNPmat.width/BSx);


fprintf(stderr,"Copying matrices to GPU success \n");
cudaEventRecord(start, 0);

GLIDE<<< dimGrid, dimBlock >>>(d_SNPmat.elements, d_SNPmat2.elements, d_Pheno.elements, d_output.elements, NSNPS_GPU, NSUB);

cudaEventRecord(stop,0);
cudaEventSynchronize(stop);
cudaEventElapsedTime(&elapsedTime,start,stop);
fprintf(stderr,"Time to do %d SNPS =%4.2f ms \n", NSNPS_GPU, elapsedTime);
fprintf(stderr,"i=%d  j=%d %d\n", index1count/NSNPS_GPU, index2count/NSNPS_GPU,NSNPS);

cudaMemcpy(h_output.elements, d_output.elements, d_mem_size_output, cudaMemcpyDeviceToHost);
fprintf(stderr,"Output Copy Success \n");


//Writing out the results
//Initialize writing file
FILE *fpout;
fpout = fopen(results, "a");

 for (int rowcounter=0;rowcounter<(d_SNPmat.width*d_SNPmat.width*NCOL_OUTPUT_GPU);rowcounter+=NCOL_OUTPUT_GPU)
  {

if(abs(h_output.elements[rowcounter+7])>=t_thres) //t-score threshold selection criterion
	{
	fprintf(fpout,"%d %d %3.0f %3.0f %3.0f %3.0f %5.5f %5.5f %5.5f %5.5f \n",
	index1count/NSNPS_GPU,
	index2count/NSNPS_GPU,
	h_output.elements[rowcounter],
	h_output.elements[rowcounter+1],
	h_output.elements[rowcounter+2],
	h_output.elements[rowcounter+3],
	h_output.elements[rowcounter+4],
	h_output.elements[rowcounter+5],
	h_output.elements[rowcounter+6],
	h_output.elements[rowcounter+7]);
	}
 }
offset2 += NSNPS_GPU*NSUB;
fclose(fpout);
cudaFreeHost(h_output.elements);

}
offset1 += NSNPS_GPU*NSUB;
offset2 = h_SNPmat2.elements; //reset 

} 



cudaFree(d_SNPmat.elements);
cudaFree(d_SNPmat2.elements);
cudaFree(d_Pheno.elements);
cudaFree(d_output.elements);

free(h_SNPmat.elements);
free(h_SNPmat2.elements);
free(h_Pheno.elements);


}

//Functions specifications
void zeroInitfloat(float* data, int size)
{
    for (int i = 0; i < size; ++i)
        data[i] = 0;
}
void randomInitInt(int* data, int size)
{
    for (int i = 0; i < size; ++i)
        data[i] = rand() / (int)RAND_MAX;
}

void parseInput(int argc, char **argv){
  int i=1;
  if(argc <= 1){
    printf("\nusage: \n  GLIDE -f1 genoFile -f2 genoFile2 -fp phenoFile -n NSubj -m NSNPS -m2 NSNPS -p NSNPS_GPU -t t_thres -o results -g deviceOrdinal \n\n");
	printf("\tgenoFile      = txt file containing the first genotype file\n");
    printf("\tgenoFile2     = txt file containing the second genotype file\n");
    printf("\tphenoFile     = txt file containing the phenotype\n");
    printf("\tNSUB			= Number of Subjects\n");
    printf("\tNSNPS			= Total Number of SNPs from genoFile\n");
    printf("\tNSNPS			= Total Number of SNPs from genoFile2\n");
	printf("\tNSNPS_GPU		= Size of the partition SNPs\n");
	printf("\tt_thres		= t-score threshold\n");
    printf("\tresults       = output file ; stored in text format\n");
    printf("\tdeviceNum     = GPU ID # use for the run\n");
	printf("\n\n");
    exit(0);
  }

  while(i<argc){
    if(!strcmp(argv[i], "-f1"))
      genoFile = argv[++i];
    else if(!strcmp(argv[i], "-f2"))
      genoFile2 = argv[++i];
    else if(!strcmp(argv[i], "-fp"))
      phenoFile = argv[++i];
    else if(!strcmp(argv[i], "-n"))
      NSUB = atoi(argv[++i]);
    else if(!strcmp(argv[i], "-m"))
      NSNPS = atoi(argv[++i]);
    else if(!strcmp(argv[i], "-m2"))
      NSNPS2 = atoi(argv[++i]);
    else if(!strcmp(argv[i], "-p"))
      NSNPS_GPU = atoi(argv[++i]);
    else if(!strcmp(argv[i], "-t"))
      t_thres = atof(argv[++i]);
	else if(!strcmp(argv[i], "-o"))
      results = argv[++i];
    else if(!strcmp(argv[i], "-g"))
      deviceOrdinal = atoi(argv[++i]);
    else{
      fprintf(stderr,"%s : unrecognized option.. exiting\n",argv[i]);
      exit(1);
    }
    i++;
  }

  if( !genoFile || !genoFile2 || !phenoFile || !NSUB || !NSNPS || !NSNPS_GPU || !results){
    fprintf(stderr,"more arguments needed.. exiting\n");
    exit(1);
  }

}

////////////////////////////////////////////////////////////////////////////////
//! Storing Matrices
//! 
////////////////////////////////////////////////////////////////////////////////
/*Storing Matrices from*/
void readfloatData(char *dataFile, unsigned int rows, unsigned int cols, float * data){
  FILE *fp;
  float *dp = data;
  int i;

  fp = fopen(dataFile,"r");
  if(fp==NULL){
    fprintf(stderr,"error opening file.. exiting\n");
    //exit(1);
  } 
  
  for (i=0; i<rows*cols; ++i){
	  fscanf(fp, "%f", dp);
	  dp++;
  } 
  fclose(fp);
}

void readintData(char *dataFile, unsigned int rows, unsigned int cols, int *data){
  FILE *fp;
  int *dp = data;
  int i;
  fp = fopen(dataFile,"r");
  if(fp==NULL){
    fprintf(stderr,"error opening file.. exiting\n");
   // exit(1);
  } 
      for (i=0; i<rows*cols; ++i){
	  fscanf(fp, "%d", dp++);
  }
  fclose(fp);
}

void readstringData(char *dataFile, unsigned int rows, unsigned int cols, char *data){
  FILE *fp;
  char *dp = data;
  int i;
  fp = fopen(dataFile,"r");
  if(fp==NULL){
    fprintf(stderr,"error opening file.. exiting\n");
   // exit(1);
  } 
      for (i=0; i<rows*cols; ++i){
	  fscanf(fp, "%s", dp++);
  }
  fclose(fp);
}

int PrintDevices(int deviceCount, int deviceSelected)
{
    cudaError_t err = cudaSuccess;

    cudaDeviceProp deviceProperty;
    for (int currentDeviceId = 0; currentDeviceId < deviceCount; ++currentDeviceId)
    {
        memset(&deviceProperty, 0, sizeof(cudaDeviceProp));
        err = cudaGetDeviceProperties(&deviceProperty, currentDeviceId);
        CheckConditionXR_(err == cudaSuccess, err);

        printf("\ndevice name: %s", deviceProperty.name);
        if (currentDeviceId == deviceSelected)
        {
            printf("    <----- creating CUcontext on this");    
        }
        printf("\n");

        printf("device sharedMemPerBlock: %Iu \n", deviceProperty.sharedMemPerBlock);
        printf("device totalGlobalMem: %Iu \n", deviceProperty.totalGlobalMem);
        printf("device regsPerBlock: %d \n", deviceProperty.regsPerBlock);
        printf("device warpSize: %d \n", deviceProperty.warpSize);
        printf("device memPitch: %Iu \n", deviceProperty.memPitch);
        printf("device maxThreadsPerBlock: %d \n", deviceProperty.maxThreadsPerBlock);
        printf("device maxThreadsDim[0]: %d \n", deviceProperty.maxThreadsDim[0]);
        printf("device maxThreadsDim[1]: %d \n", deviceProperty.maxThreadsDim[1]);
        printf("device maxThreadsDim[2]: %d \n", deviceProperty.maxThreadsDim[2]);
        printf("device maxGridSize[0]: %d \n", deviceProperty.maxGridSize[0]);
        printf("device maxGridSize[1]: %d \n", deviceProperty.maxGridSize[1]);
        printf("device maxGridSize[2]: %d \n", deviceProperty.maxGridSize[2]);
        printf("device totalConstMem: %Iu \n", deviceProperty.totalConstMem);
        printf("device major: %d \n", deviceProperty.major);
        printf("device minor: %d \n", deviceProperty.minor);
        printf("device clockRate: %d \n", deviceProperty.clockRate);
        printf("device textureAlignment: %Iu \n", deviceProperty.textureAlignment);
        printf("device deviceOverlap: %d \n", deviceProperty.deviceOverlap);
        printf("device multiProcessorCount: %d \n", deviceProperty.multiProcessorCount);

        printf("\n");
    }

    return cudaSuccess;
}
