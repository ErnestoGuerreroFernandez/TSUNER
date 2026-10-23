.SUFFIXES:	.cu .c

export CUDA_PATH=/usr/local/cuda
export UNAME=$(shell uname -m)

# Compilers
export NVCC 		= $(CUDA_PATH)/bin/nvcc
export CXX 		= g++-14

DEPS = ./GPU/libkernels.a

export CXXFLAGS         =-O2
export CXXINC          =-I./ -I$(CUDA_PATH)/include -I/share/apps/NETCDF_C/include
export CXXLIBS         =-L./GPU  -lkernels -L$(CUDA_PATH)/lib64 -L/usr/local/cuda/lib64 -lcudart -L/share/apps/NETCDF_C/lib -lnetcdf
OBJSC 	:= main_TSUNER.o

%.o: 	%.cxx
	$(CXX) $(CXXFLAGS) $(CXXINC) -c $*.cxx


all: kernels tsuner.exe
kernels:
	cd GPU; make; cd ..;

tsuner.exe: $(OBJSC)
	$(CXX)  -o $@ $(OBJSC) $(CXXLIBS)  

.PHONY: clean
clean:
	rm -fr *.o *~ *.exe
	cd GPU; make clean; cd ..;
