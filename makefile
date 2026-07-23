# Compilers
CC       = gcc
CXX      = g++
NVCC     = nvcc
TARGET   = exe.exe
BUILDDIR = build
SRCDIR   = src

# Change the below paths so it corresponds to your system
MAT_NUM_INC = /home/erdi/Uni/MyPrograms/MatLib/MathNum/Include
LAGJET_INC_DIR_SRC = /home/erdi/Uni/LCGAStage/Jets/LagrangianJet/Include

# Flags
INCLUDES = -I. -I$(MAT_NUM_INC) -I$(LAGJET_INC_DIR_SRC)
# Base flags for both languages
BASE_FLAGS = -pg -Wall -Wextra -Wshadow -Wcast-align -Wfloat-equal -Wundef -Wcast-qual $(INCLUDES) -O3 -fopenmp
CFLAGS     = $(BASE_FLAGS) 
CXXFLAGS   = $(BASE_FLAGS) -fpermissive
CUFLAGS    = -O3 $(INCLUDES)
LDFLAGS    = -lm -lstdc++ -fopenmp

# Sources and Objects
# Find all .c and .cpp in SRCDIR and MAT_NUM_INC
SRC_C   = $(wildcard $(SRCDIR)/*.c)   $(shell find $(MAT_NUM_INC) -name '*.c')
SRC_CPP = $(wildcard $(SRCDIR)/*.cpp) $(wildcard $(LAGJET_INC_DIR_SRC)/*.cpp)
SRC_CU  = $(wildcard $(SRCDIR)/*.cu)

# Map sources to object files in the BUILDDIR
OBJS = $(addprefix $(BUILDDIR)/, $(notdir $(SRC_C:.c=.o))) \
       $(addprefix $(BUILDDIR)/, $(notdir $(SRC_CPP:.cpp=.o))) \
       $(addprefix $(BUILDDIR)/, $(notdir $(SRC_CU:.cu=.o)))

# VPATH allows Make to find the source files regardless of directory
VPATH = $(SRCDIR):$(shell find $(MAT_NUM_INC) -type d):$(shell find $(LAGJET_INC_DIR_SRC) -type d)

all: $(BUILDDIR) $(TARGET)

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

# Link step: Use nvcc to link if there is any cuda code
$(TARGET): $(OBJS)
	$(CXX) $(BASE_FLAGS) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)


# $(NVCC) $(CUFLAGS) -o $@ $^ $(LDFLAGS)

# Rule for C files
$(BUILDDIR)/%.o: %.c
	@echo "Compiling C file: $<..."
	$(CC) $(CFLAGS) -c $< -o $@ $(LDFLAGS)

# Rule for C++ files
$(BUILDDIR)/%.o: %.cpp
	@echo "Compiling C++ file: $<..."
	$(CXX) $(CXXFLAGS) -c $< -o $@ $(LDFLAGS)


# Rule for CUDA files
$(BUILDDIR)/%.o: %.cu
	@echo "Compiling C++ file: $<..."
	$(NVCC) $(CUFLAGS) -c $< -o $@


# Utility targets
.PHONY: clean debug run

debug: CFLAGS += -g -O0 -fsanitize=address
debug: CXXFLAGS += -g -O0 -fsanitize=address
debug: LDFLAGS += -fsanitize=address
debug: clean all

run: all
	./$(TARGET)

clean:
	rm -rf $(BUILDDIR) $(TARGET)