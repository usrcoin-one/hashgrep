CC=nvcc
CXX=nvcc
OPT=-O3
#OPT=-g -ggdb
CFLAGS=-c $(OPT)
NVCFLAGS=-gencode arch=compute_20,code=sm_20 
NVLDFLAGS= #-L./SDKlink/lib -lcutil_i386
HEADERS = $(wildcard *.h)
TARGETS = greptest hashgrep

all: $(TARGETS)

greptest:  greptest.o getfile.o timing.o
	nvcc -o greptest -gencode arch=compute_20,code=sm_20 greptest.o getfile.o timing.o $(LDFLAGS)

hashgrep:  hashgrep.o getfile.o timing.o
	nvcc -o hashgrep -gencode arch=compute_20,code=sm_20 hashgrep.o getfile.o timing.o $(LDFLAGS)


greptest.o: greptest.cu $(HEADERS) Makefile
	nvcc $(NVCFLAGS) -c -o greptest.o greptest.cu -I./SDKlink/common/include -g

hashgrep.o: hashgrep.cu $(HEADERS) Makefile
	nvcc $(NVCFLAGS) -c -o hashgrep.o hashgrep.cu -I./SDKlink/common/include -g

%.o: %.c
	gcc $(CFLAGS) $< -o $@  -g

clean:
	/bin/rm -f *.o *~ greptest hashgrep
