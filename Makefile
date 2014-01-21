VPATH=src/
OBJ=cgaus.o cgmat.o cgif.o
CC=gcc
CFLAGS=-std=c99
LFLAGS=

all: cgaus

cgaus: ${OBJ}
	${CC} ${CFLAGS} ${LFLAGS} -o cgaus $+

${OBJ}: %.o: %.c
	${CC} -c ${CFLAGS} ${LFLAGS} -o $@ $<

${OBJ}: cgaus.h cgmat.h

clean:
	rm -f ${OBJ} cgaus

.PHONY: all clean	
