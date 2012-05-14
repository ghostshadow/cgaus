VPATH=src/
OBJ=cgaus.o cgmat.o cgif.o
CC=cc
CFLAGS=-std=c99

all: cgaus

cgaus: -lncurses ${OBJ}
	${CC} ${CFLAGS} -o cgaus $+

${OBJ}: %.o: %.c
	${CC} -c ${CFLAGS} -o $@ $<

${OBJ}: cgaus.h

clean:
	rm -f ${OBJ}

.PHONY: all clean	
