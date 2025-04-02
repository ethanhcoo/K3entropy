
CFLAGS = -g -I/opt/homebrew/include #changed from '-g'
LDFLAGS = -L/opt/homebrew/lib -lmpfr -lgmp #-lmpfr -lgmp used for precision calculations
BIN    = ../bin
CC     = gcc

all:	stable mclass periodic periodicExact derivatives derivativesExact 

stable:	stable.o k3.o ps.o 
	${CC} stable.o k3.o ps.o -lm -o ${BIN}/stable

mclass: mclass.o k3.o ps.o 
	${CC} mclass.o k3.o ps.o -lm -o ${BIN}/mclass

periodic: periodic.o k3.o ps.o 
	${CC} periodic.o k3.o ps.o ${LDFLAGS} -lm -o ${BIN}/periodic

periodicExact: periodicExact.o k3.o ps.o 
	${CC} periodicExact.o k3.o ps.o ${LDFLAGS} -lm -o ${BIN}/periodicExact
	
derivatives: derivatives.o k3.o ps.o 
	${CC} derivatives.o k3.o ps.o ${LDFLAGS} -lm -o ${BIN}/derivatives

derivativesExact: derivativesExact.o k3.o ps.o 
	${CC} derivativesExact.o k3.o ps.o ${LDFLAGS} -lm -o ${BIN}/derivativesExact

clean:
	rm -f *.o ${BIN}/*
