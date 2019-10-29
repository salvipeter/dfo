all: test

CXXFLAGS=-Wall -pedantic -std=c++17 -O3

test: test.o downhill.o powell.o
	g++ -o $@ $^
