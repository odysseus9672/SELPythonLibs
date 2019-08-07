SOURCES=CosmoCalcs.c test.c
OBJECTS=$(SOURCES:.c=.o)
SHARED=libCosmoCalcs.so
EXECUTABLES=test

CFLAGS=$(CFLAGS0)
LDFLAGS=$(LDFLAGS0) -lgsl -lgslcblas

CC=gcc

all: shared $(OBJECTS) $(EXECUTABLES)

shared: CosmoCalcs.o
	$(CC) $(CFLAGS) -shared -o libCosmoCalcs.so CosmoCalcs.o $(LDFLAGS)

clean:
	rm $(OBJECTS) $(SHARED) $(EXECUTABLES)

CosmoCalcs.o:

test.o:

test: test.o CosmoCalcs.o
	$(CC) $(CFLAGS) $(OBJECTS) -o test $(LDFLAGS)