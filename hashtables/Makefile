
UNAME := $(shell uname -a)
ifeq ($(firstword $(filter Linux,$(UNAME))),Linux)
  CC = g++
endif
ifeq ($(firstword $(filter Darwin,$(UNAME))),Darwin)
  CC = g++-mp-4.5
endif

# Uncomment one of the following to switch between debug and opt mode
#OPT = -std=c++0x -O3 -DNDEBUG
OPT = -std=c++0x -g -ggdb

CFLAGS += -Wall -c -I. -I./include -I/usr/include/ $(OPT)

LDFLAGS+= -Wall 
#-lpthread -lssl -lcrypto

LIBOBJECTS = \
#	./CSA.o \


HEADERS = $(wildcard *.h)

TESTS = \
	hashtable_test\

BENCH = \

#	filter_bench \
	mem_bench


PROGRAMS = $(TESTS) $(BENCH)

all: $(PROGRAMS)

check: $(TESTS)
	for t in $(TESTS); do echo "***** Running $$t"; ./$$t || exit 1; done

clean:
	rm -f $(PROGRAMS) *.o */*.o

cleanall: clean
	rm -f */*.gcda

hashtable_test: hashtable_test.o $(LIBOBJECTS) 
	$(CC) hashtable_test.o $(LIBOBJECTS) $(LDFLAGS) -o $@


%.o: %.cc ${HEADERS} Makefile
	$(CC) $(CFLAGS) $< -o $@

# TODO(gabor): dependencies for .o files
# TODO(gabor): Build library
