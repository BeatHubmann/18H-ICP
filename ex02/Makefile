CXX = g++-8
CXXFLAGS = -std=c++11 -O3 -march=native -funroll-loops -fopenmp

.PHONY: all
all: example

example: example.cpp

task1: task1.cpp

task2: task2.cpp

image2: task2
	./task2 0.7 100
	convert task2_*.png task2.gif

.PHONY: clean
clean:
	@rm -rf example task1 task2 task3 *.ppm *.png *.gif