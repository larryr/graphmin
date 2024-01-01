CC = clang

help:
	@echo all -- build directory into a static library
	@echo clean -- remove generated files

# Get all .c files in this directory
SRCS = $(wildcard *.c)

# Define object files
OBJS = $(SRCS:.c=.o) gen-grammar.o gen-scan.o

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

# Rule to build .c file from .y file.
gen-%.c gen-%.h: %.y
	bison --defines=gen-$*.h --output=gen-$*.c $<

# Rule to build .c .h files from .l file.
gen-%.c gen-%.h: %.l
	flex --outfile=gen-$*.c --header-file=gen-$*.h $<


.PHONY: all clean

clean:
	rm -f $(OBJS) $(LIB) gen-*
	z