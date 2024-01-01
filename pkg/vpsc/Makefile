CC = clang

help:
	@echo all -- build directory into a static library
	@echo clean -- remove generated files

# Get all .c files in this directory
SRCS = $(wildcard *.cpp)

# Define object files
OBJS = $(SRCS:.cpp=.o)

# Define static library name
DIR_NAME = $(notdir $(CURDIR))
LIB = lib$(DIR_NAME).a

# Defining compiler flags
CFLAGS = -I.. -I. -std=c++20 -stdlib=libc++

# Default rule
all: $(LIB)

# Rule to build the .a file.
$(LIB): $(OBJS)
	ar rcs $@ $^

# Rules to build .o file from .cpp file.
%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: all clean

clean:
	rm -f $(OBJS) $(LIB)