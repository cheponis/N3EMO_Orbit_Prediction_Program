
# 1. Define variables
CC      = cc
CFLAGS  = -g -Wall -Wextra -std=gnu2x -D_NETBSD_SOURCE -Wformat=2 \
          -Wswitch-enum -Wshadow -Wstrict-prototypes -Wold-style-definition \
          -Wvla -Wc11-c2x-compat -O2
LIBS    = -lm

# 2. Define targets
all: orbit

# 3. Linking the executable
orbit: orbit.o orbitr.o
        $(CC) $(CFLAGS) orbit.o orbitr.o -o orbit $(LIBS)

# 4. Compiling source files into object files
orbit.o: orbit.c
        $(CC) $(CFLAGS) -c orbit.c

orbitr.o: orbitr.c
        $(CC) $(CFLAGS) -c orbitr.c

# 5. Clean up utility
.PHONY: clean
clean:
        rm -f *.o orbit
