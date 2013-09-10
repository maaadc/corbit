CC=g++
CFLAGS=-Wall -O2
LDFLAGS=

TARGET=corbit

SRCDIR=src
BUILDDIR=build

SOURCES=$(shell find $(SRCDIR) -type f -name *.cc)
OBJECTS=$(SOURCES:.cpp=.o)

all: $(SOURCES) $(TARGET)

$(TARGET): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

$(SOURCES):
	$(CC) $(CFLAGS) $< -o $@

clean:
	$(RM) -r $(BUILDDIR) $(TARGET)
