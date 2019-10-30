all: test

CXXFLAGS=-Wall -pedantic -std=c++17 -O3

test: test.o direct.o mads.o nelder-mead.o powell.o
	g++ -o $@ $^
