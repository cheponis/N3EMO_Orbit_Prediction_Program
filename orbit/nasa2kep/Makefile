# Define compiler and flags
CC      = cc
CFLAGS  = -g -Wall -Wextra -std=gnu2x -Wformat=2 \
          -Wswitch-enum -Wshadow -Wstrict-prototypes -Wold-style-definition \
          -Wvla -Wc11-c2x-compat -O2
LIBS    = -lm

# Default target
all: nasa2kep

# Linking the final executable
# Note: $(LIBS) comes at the end
nasa2kep: nasa2kep.o
	$(CC) $(CFLAGS) nasa2kep.o -o nasa2kep $(LIBS)

# Compiling the object file
nasa2kep.o: nasa2kep.c
	$(CC) $(CFLAGS) -c nasa2kep.c

# Clean up
.PHONY: clean
clean:
	rm -f *.o nasa2kep
	[ -f kepler.dat ] && mv kepler.dat kepler.dat.bak || true
