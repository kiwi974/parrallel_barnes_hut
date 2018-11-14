#OPTIM+=-O3 -march=native
DEBUG+=-g -O1
CXX=g++
CC=g++
LD=${CXX}
CXXFLAGS+=$(OPTIM) $(DEBUG) -Wall -Wextra -std=c++11 -pg -O1 -no-pie
LDFLAGS+=$(OPTIM) $(DEBUG) $(CXXFLAGS) -lm -pg -no-pie

all: main

pi: main.o

clean:
	rm -f main *.o *~
