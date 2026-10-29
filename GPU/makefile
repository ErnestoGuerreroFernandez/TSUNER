.SUFFIXES: .cu

export CUDA_PATH ?= /usr/local/cuda
ARCH ?= native


export NVCC      = $(CUDA_PATH)/bin/nvcc
export INC       = -I./ -I../ -I/share/apps/NETCDF_C/include/
export LIB       = -L$(CUDA_PATH)/lib64
export NVCCFLAGS = -O2 -arch=$(ARCH) -fmad=false -rdc=true -dlto -MMD -MP

# Directory for auxiliary GPU build files
OBJ_DIR := obj

SRCS := CheckCUDASupport.cu InitializeGPU.cu netcdf.cu Reconstructor.cu Reduction_kernel_max.cu Reduction_kernel_min.cu main_TSUNER.cu Volume_kernel.cu VolumeSolver.cu

# Place all GPU .o files inside obj/
OBJS := $(patsubst %.cu, $(OBJ_DIR)/%.o, $(SRCS))

all: libkernels.a

# Compile .cu into obj/%.o
$(OBJ_DIR)/%.o: %.cu | $(OBJ_DIR)
	$(NVCC) $(INC) $(LIB) $(NVCCFLAGS) -c $< -o $@

# Create obj/ folder if it doesn't exist
$(OBJ_DIR):
	mkdir -p $@

# Intermediate device linking placed cleanly inside obj/
$(OBJ_DIR)/dlink.o: $(OBJS)
	$(NVCC) $(LIB) $(NVCCFLAGS) -dlink $(OBJS) -o $@

# Pack static library
libkernels.a: $(OBJS) $(OBJ_DIR)/dlink.o
	$(NVCC) --lib $(OBJS) $(OBJ_DIR)/dlink.o -o libkernels.a

# Include dependency files from obj/
-include $(OBJS:.o=.d)

.PHONY: clean
clean:
	rm -rf $(OBJ_DIR) *.o *.d *~ libkernels.a
