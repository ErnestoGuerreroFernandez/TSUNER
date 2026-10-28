.SUFFIXES: .cu .c .cxx

export CUDA_PATH = /usr/local/cuda
export UNAME     = $(shell uname -m)

# Compilers
export NVCC      = $(CUDA_PATH)/bin/nvcc
export CXX       = g++

# Directory for auxiliary build files
OBJ_DIR := obj

# Source files
SRCS    := main_TSUNER.cxx Problem_loader.cxx interpolations.cxx cini.cxx tides.cxx

# Place all .o files inside $(OBJ_DIR)/
OBJSC   := $(patsubst %.cxx, $(OBJ_DIR)/%.o, $(SRCS))

export CXXFLAGS  = -O2 -MMD -MP
export CXXINC    = -I./ -I$(CUDA_PATH)/include -I/share/apps/NETCDF_C/include
export CXXLIBS   = -L./GPU -Wl,--start-group -lkernels -Wl,--end-group \
                   -L$(CUDA_PATH)/lib64 -lcudart \
                   -L/share/apps/NETCDF_C/lib -lnetcdf

all: kernels tsuner.exe

# Compile .cxx into obj/%.o (the '|' ensures the directory is created first)
$(OBJ_DIR)/%.o: %.cxx | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) $(CXXINC) -c $< -o $@

# Create obj/ folder if it doesn't exist
$(OBJ_DIR):
	mkdir -p $@

# Build the GPU library inside GPU/
kernels:
	$(MAKE) -C GPU

# Link final executable
tsuner.exe: $(OBJSC) GPU/libkernels.a
	$(CXX) -o $@ $(OBJSC) $(CXXLIBS)

# Include dependency files located in obj/
-include $(OBJSC:.o=.d)

.PHONY: clean kernels
clean:
	rm -rf $(OBJ_DIR) *.o *.d *~ *.exe
	$(MAKE) -C GPU clean
