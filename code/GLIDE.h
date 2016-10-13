/* This file is part of the EpiGPULin (EGL) library.
* (C) Copyright 2011, Tony Kam-Thong [tony@mpipsykl.mpg.de]
*/

#ifndef LinReg_H
#define LinReg_H

/* the BSx values are replaced by shell script */

typedef struct {
int width;
int height;
int* elements;
} Matrix;

typedef struct {
int width;
int height;
float* elements;
} Matrixfloat;

//Setting BLOCK_SIZE
/*
//sm_20
#define BSx 16
#define BSy 16
#define BSz 1
//sm_20
*/

//begin:sm_13
#define BSx 10
#define BSy 10
#define BSz 1
//end:sm_13



#define p 10 //number of subjects info read in at a time, arbitrary for now;

#define IDX(i,j,ld) (((i)*(ld))+(j))

#define NPAR 5

#define NCOL_OUTPUT 10
#define NCOL_OUTPUT_GPU 8

#endif
