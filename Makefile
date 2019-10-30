all: test

CXX=g++-9
CXXFLAGS=-Wall -pedantic -std=c++17 -O3

test: test.o mads.o nelder-mead.o powell.o
	g++ -o $@ $^
