CC:=gcc
CFLAGS:=-O3 -Wall -g -pg

.c.o:
	$(CC) -c $(CFLAGS) $<

.PHONY: all
all:    navier

.PHONY: clean
clean:
	rm -f *.o navier

navier: alloc.o boundary.o init.o main.o output.o simulation.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

boundary.o       : datadef.h
init.o           : datadef.h
main.o           : alloc.h boundary.h datadef.h init.h simulation.h
simulation.o     : datadef.h init.h

