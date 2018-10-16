CXX = g++-8
CXXFLAGS = -std=c++1z -O3 -march=native -funroll-loops 

.PHONY: all
all: task1 task2 figure1 figure2 

task1: task1.cpp

figure1: task1 figure1.gplt
	./task 10 1000 1 100 > results.dat
	gnuplot figure.gplt

.PHONY: clean
clean:
	@rm -rf task *.png *.dat 