
CFLAGS = -g -I/opt/homebrew/include #changed from -g
LDFLAGS = -L/opt/homebrew/lib -lgmp #added this
BIN    = ../bin
CC     = gcc

all:	stable mclass periodic derivatives# traj

stable:	stable.o k3.o ps.o 
	${CC} stable.o k3.o ps.o -lm -o ${BIN}/stable

mclass: mclass.o k3.o ps.o 
	${CC} mclass.o k3.o ps.o -lm -o ${BIN}/mclass

periodic: periodic.o k3.o ps.o 
	${CC} periodic.o k3.o ps.o ${LDFLAGS} -lm -o ${BIN}/periodic
	
derivatives: derivatives.o k3.o ps.o 
	${CC} derivatives.o k3.o ps.o -lm -o ${BIN}/derivatives

#traj:   traj.o k3.o cx.o scan.o 
#	${CC} traj.o k3.o cx.o scan.o -lm -o ${BIN}/traj

