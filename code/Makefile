# This file is part of the GLIDE library.
# (C) Copyright 2011, Tony Kam-Thong [tony@mpipsykl.mpg.de]
CC=gcc
NVCC=nvcc
CCFLAGS=
# NVCCFLAGS= -arch=sm_20

# Older CUDA GPU with compute capability of 1.3
NVCCFLAGS= -arch=sm_13  --compiler-bindir=/usr/bin/gcc-4.3

# other flags: -deviceemu -arch=sm_20 --ptxas-options=-v
# These are useful when debugging sometimes.
LINKFLAGS=-lcuda  -lm 
# other linkflags: -lgsl -lgslcblas -lcublas
# lgsl and lgslcblas are required if you want to use the GSL.
# SOURCES=
CUSOURCES=GLIDE_main.cu
OBJECTS=$(SOURCES:.c=.o)
CUOBJECTS=$(CUSOURCES:.cu=.o)
EXECUTABLE= GLIDE
all: $(SOURCES) $(CUSOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) $(CUOBJECTS)
	$(NVCC) $(NVCCFLAGS) $(OBJECTS) $(CUOBJECTS) -o $@ $(LINKFLAGS)

%.o:%.c
	$(NVCC) $(NVCCFLAGS) -c $+

%.o:%.cu
	$(NVCC) $(NVCCFLAGS) -c $+

clean:
	rm -f *.o
	rm -f $(EXECUTABLE)

GLIDE_main.cu: GLIDE_kernel_op.cu
	echo kernel changed

# END
