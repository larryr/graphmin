
CC = clang
CLANG = clang++

# Get all .c and .cpp files in this directory
SRCS = $(wildcard *.c *.cpp)

# Define object files for both .c and .cpp files
OBJS = $(patsubst %.c,%.o,$(patsubst %.cpp,%.o,$(SRCS)))

# Get current directory name
DIR_NAME = $(notdir $(CURDIR))

# Library base names
LIB_BASE_NAMES = cdt cgraph circogen common dotgen fdpgen gvc label neatogen \
ortho osage pack patchwork pathplan rbtree sparse twopigen xdot vpsc sfdpgen \
plugin_core plugin_dot_layout plugin_neato_layout

LIB_SYS = expat z

# Construct library search paths for the linker (-L)
LIB_PATHS = $(addprefix -L../../pkg/,$(LIB_BASE_NAMES))

# Construct library names for the linker (-l)
LIB_ALL_NAMES = $(LIB_BASE_NAMES) $(LIB_SYS)
LIB_NAMES = $(addprefix -l,$(LIB_ALL_NAMES))

# Define executable name as current directory
EXEC = $(DIR_NAME)

# Defining compiler flags
CFLAGS = -I. -I../../pkg

# Default rule
all: $(EXEC)

# Rule to build the executable
$(EXEC): $(OBJS)
	$(CLANG) -v -o $@ $^ $(LIB_PATHS) $(LIB_NAMES)

# Rules to build .o file from .c file.
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Rule to build .o file from .cpp file.
%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: all clean

clean:
	rm -f $(OBJS) $(EXEC)