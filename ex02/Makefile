CXX = g++-8
CXXFLAGS = -std=c++1z -O3 -march=native -funroll-loops -fopenmp

.PHONY: all
all: example task1 image2 task3

example: example.cpp

task1: task1.cpp

task2: task2.cpp

image2: task2
	./task2 0.6 100 10 1
	convert task2_*.png task2.gif

task3: task3.cpp

.PHONY: clean
clean:
	@rm -rf example task1 task2 task3 *.ppm *.png *.gif