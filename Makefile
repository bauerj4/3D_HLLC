CONFIG   =  Config.sh
PYTHON   =  /usr/bin/python
CFLAGS   =
CC       = gcc
INCL	 = include/allvars.h include/hllc_defs.h include/proto.h
SRC 	 = src/main.c src/mesh.c src/flux.c src/eos.c src/bcs.c src/formulation.c

make:
	python make_macros.py
	$(CC) $(CFLAGS) $(SRC) $(INCL) -Wall -O3 -lm -o bin/hllc3

clean:
	rm *~ src/*.o src/*~ include/*.o include/*~ bin/*
