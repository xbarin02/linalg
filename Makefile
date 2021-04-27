CXX=g++
CXXFLAGS=-std=c++11 -pedantic -Wall -Wextra -Og -g
LDFLAGS=-rdynamic
LINK.o=g++
BIN=main

.PHONY: all
all: $(BIN)

main.o: main.cpp Vector.h Matrix.h

.PHONY: clean
clean:
	$(RM) -- *.o $(BIN)
