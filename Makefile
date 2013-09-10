CC=g++
CFLAGS=-Wall -O2
LDFLAGS=
SOURCES=corbit.cc data.cc sim.cc
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=corbit

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@
