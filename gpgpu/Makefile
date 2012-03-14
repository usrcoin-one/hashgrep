CC=nvcc
CXX=nvcc
LDFLAGS= #-L./SDKlink/lib -lcutil_i386

all: greptest hashgrep 

greptest:  greptest.o getfile.o timing.o
	nvcc -o greptest -gencode arch=compute_20,code=sm_20 greptest.o getfile.o timing.o $(LDFLAGS)

hashgrep:  hashgrep.o getfile.o timing.o
	nvcc -o hashgrep -gencode arch=compute_20,code=sm_20 hashgrep.o getfile.o timing.o $(LDFLAGS)


greptest.o: greptest.cu
	nvcc -gencode arch=compute_20,code=sm_20 -c -o greptest.o greptest.cu -I./SDKlink/common/include -g

hashgrep.o: hashgrep.cu
	nvcc -gencode arch=compute_20,code=sm_20 -c -o hashgrep.o hashgrep.cu -I./SDKlink/common/include -g


getfile.o: getfile.c
	gcc $(CFLAGS) -c getfile.c -g

clean:
	/bin/rm -f *.o *~ greptest hashgrep
