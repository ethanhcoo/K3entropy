
CFLAGS = -g -I/opt/homebrew/include #changed from '-g'
LDFLAGS = -L/opt/homebrew/lib -lmpfr -lgmp #-lmpfr -lgmp used for precision calculations
BIN    = ../bin
CC     = gcc

all:	stable mclass periodic RecurrentOrbit derivatives# traj

stable:	stable.o k3.o ps.o 
	${CC} stable.o k3.o ps.o -lm -o ${BIN}/stable

mclass: mclass.o k3.o ps.o 
	${CC} mclass.o k3.o ps.o -lm -o ${BIN}/mclass

periodic: periodic.o k3.o ps.o 
	${CC} periodic.o k3.o ps.o ${LDFLAGS} -lm -o ${BIN}/periodic

RecurrentOrbit: RecurrentOrbit.o k3.o ps.o 
	${CC} RecurrentOrbit.o k3.o ps.o ${LDFLAGS} -lm -o ${BIN}/RecurrentOrbit
	
derivatives: derivatives.o k3.o ps.o 
	${CC} derivatives.o k3.o ps.o ${LDFLAGS} -lm -o ${BIN}/derivatives

clean:
	rm -f *.o ${BIN}/*
#traj:   traj.o k3.o cx.o scan.o 
#	${CC} traj.o k3.o cx.o scan.o -lm -o ${BIN}/traj

