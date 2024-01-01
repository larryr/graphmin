CC = clang

help:
	@echo all -- build directory into a static library
	@echo clean -- remove generated files

# Get all .c files in this directory
SRCS = $(wildcard *.c)

# Define object files
OBJS = $(SRCS:.c=.o)

# Define static library name
DIR_NAME = $(notdir $(CURDIR))
LIB = lib$(DIR_NAME).a

# Defining compiler flags
CFLAGS = -I.. -I.

# Default rule
all: $(LIB)

# Rule to build the .a file.
$(LIB): $(OBJS)
	ar rcs $@ $^

# Rules to build .o file from .c file.
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY: all clean

clean:
	rm -f $(OBJS) $(LIB)