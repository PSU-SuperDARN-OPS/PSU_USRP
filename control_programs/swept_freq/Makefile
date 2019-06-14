OS:=$(shell uname -s)

CCpp=g++
CFLAGS=-g
LFLAGS=-L"/usr/lib/x86_64-linux-gnu/hdf5/serial"
libs=-lboost_program_options -lhdf5
SOURCES=main.cpp
OBJECTS=$(SOURCES:.cpp=.o)
#EXECUTABLE=${SOUNDER_DIR}/build/swept_freq
EXECUTABLE=swept_freq

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CCpp) $(CFLAGS) -o $@ $(OBJECTS) $(libs) $(LFLAGS)

.cpp.o:
	$(CCpp) $(CFLAGS) $< -o $@ -c

clean:
	rm -rf *.o client

