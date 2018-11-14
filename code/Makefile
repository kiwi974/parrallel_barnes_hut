#OPTIM+=-O3 -march=native
DEBUG+=-g -O1
CXX=mpicxx
CC=mpicxx
LD=${CXX}
CXXFLAGS+=$(OPTIM) $(DEBUG) -Wall -Wextra -std=c++11 -pg -O1 -no-pie -fopenmp
LDFLAGS+=$(OPTIM) $(DEBUG) $(CXXFLAGS) -lm -pg -no-pie

all: mainMPI

pi: mainMPI.o

clean:
		rm -f mainMPI *.o *~
