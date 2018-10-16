CXX = g++-8
CXXFLAGS = -std=c++1z -O3 -march=native -funroll-loops 

.PHONY: all
all: task figure 

task: task.cpp

figure: task figure.gplt
	./task 10 1000 1 100 > results.dat
	gnuplot figure.gplt

.PHONY: clean
clean:
	@rm -rf task *.png *.dat 