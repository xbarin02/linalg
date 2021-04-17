CXX=g++
CXXFLAGS=-std=c++11 -pedantic -Wall -Wextra
LINK.o=g++
BIN=main

.PHONY: all
all: $(BIN)

main.o: main.cpp Vector.h

.PHONY: clean
clean:
	$(RM) -- *.o $(BIN)
