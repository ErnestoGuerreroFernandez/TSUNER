.SUFFIXES:	.cu


export CUDA_PATH=/usr/local/cuda

export NVCC	=$(CUDA_PATH)/bin/nvcc
export INC      = -I../ -I/share/apps/NETCDF_C/include/
export LIB	= -L$(CUDA_PATH)/lib64

export NVCCFLAGS= -fmad=false 

OBJS	:= Edge_kernel.o CheckCUDASupport.o InitializeGPU.o SedFlux.o netcdf.o Reconstructor.o Reduction_kernel_max.o Reduction_kernel_min.o main_TSUNER.o Volume_kernel.o VolumeSolver.o

%.o : %.cu
	$(NVCC) $(INC) $(LIB) $(NVCCFLAGS) -c $*.cu
	
all: libkernels.a

libkernels.a : $(OBJS)
	$(NVCC) --lib $(OBJS) -o libkernels.a	
	
.PHONY: clean
clean:
	rm -fr *.o *~
	rm libkernels.a
