#CC := /Developer/NVIDIA/CUDA-8.0/bin/nvcc
#CC := /usr/local/cuda/bin/nvcc
CC := nvcc
#CC := g++ # This is the main compiler
#CC := clang --analyze # and comment out the linker last line for sanity
SRCDIR := src
BUILDDIR := build
TARGETDIR := bin
TARGET := bin/runner
 
C_SRC_EXT := c
CPP_SRC_EXT := cpp
CU_SRC_EXT := cu

C_SOURCES := $(shell find $(SRCDIR) -type f -name *.$(C_SRC_EXT))
CPP_SOURCES := $(shell find $(SRCDIR) -type f -name *.$(CPP_SRC_EXT))
CU_SOURCES := $(shell find $(SRCDIR) -type f -name *.$(CU_SRC_EXT))

C_OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(C_SOURCES:.$(C_SRC_EXT)=.o))
CPP_OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(CPP_SOURCES:.$(CPP_SRC_EXT)=.o))
CU_OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(CU_SOURCES:.$(CU_SRC_EXT)=.o))

SRC_DIRS := $(shell find $(SRCDIR) -type d)

BUILD_PRE := build/
SRC_PRE := src/

BUILD_DIRS := $(subst $(SRC_PRE),$(BUILD_PRE),$(SRC_DIRS))

OBJECTS := $(C_OBJECTS) $(CPP_OBJECTS) $(CU_OBJECTS)

#-g = Host Debug...
#-G = Device Debug...

DEBUG_CFLAGS := -g -G -std=c++11 
PROD_CFLAGS := -O3 -std=c++11 

#CFLAGS := $(DEBUG_CFLAGS)
CFLAGS := $(PROD_CFLAGS)

# OpenGL+GLUT OS-specific define
ifeq ($(shell uname -s),Darwin)
GLUT_LIBS := -Xlinker -framework -Xlinker GLUT -Xlinker -framework -Xlinker OpenGL
else
GLUT_LIBS := -lGL -lGLU -lglut 
endif

CU_FLAGS := --ptxas-options=-v -gencode arch=compute_61,code=compute_61 -gencode arch=compute_61,code=sm_61
#CU_FLAGS := 

LIB := 
INC := -I include -I cuda-samples/Common// 

$(TARGET): $(OBJECTS)
	@echo 'Building target: $@'
	@echo 'Invoking: NVCC Linker'
	$(CC) --cudart static $(LIB) $(GLUT_LIBS) --relocatable-device-code=true $(CU_FLAGS) -link -o  $(TARGET) $(OBJECTS)
	@if [ $(shell uname -s) = "Darwin" ]; then dsymutil $@; fi
	@echo 'Finished building target: $@'
	@echo ' '

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(C_SRC_EXT)
	@echo 'Building C file: $<'
	@echo 'Invoking: NVCC Compiler'
	$(CC) $(INC) $(CFLAGS) $(CU_FLAGS) -M -o $@ $<
	$(CC) $(INC) $(CFLAGS) --compile  -x c++ -o  $@ $<
	@echo 'Finished building: $<'
	@echo ' '

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(CPP_SRC_EXT)
	@echo 'Building CPP file: $<'
	@echo 'Invoking: NVCC Compiler'
	$(CC) $(INC) $(CFLAGS) $(CU_FLAGS) -M -o $@ $<
	$(CC) $(INC) $(CFLAGS) --compile  -x c++ -o  $@ $<
	@echo 'Finished building: $<'
	@echo ' '

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(CU_SRC_EXT)
	@echo 'Building CU file: $<'
	@echo 'Invoking: NVCC Compiler'
	$(CC) $(INC) $(CFLAGS) $(CU_FLAGS) -M -o "$@" "$<"
	$(CC) $(INC) $(CFLAGS) --compile --relocatable-device-code=true $(CU_FLAGS)  -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

all: $(TARGET)
	@echo "Compilation Completed!"

clean:
	@echo " Cleaning..."; 
	@echo " $(RM) -r $(BUILDDIR)";
	@echo " $(RM) $(TARGET)"; 
	$(RM) -r $(BUILDDIR)
	$(RM) $(TARGET)
	@mkdir -p $(BUILDDIR)
	@mkdir -p $(BUILD_DIRS)
	@mkdir -p $(TARGETDIR)
	@echo $(BUILD_DIRS)

# Tests
tester:
	$(CC) $(CFLAGS) test/tester.cpp $(INC) $(LIB) -o bin/tester

# Spikes
ticket:
	$(CC) $(CFLAGS) spikes/ticket.cpp $(INC) $(LIB) -o bin/ticket

.PHONY: clean

