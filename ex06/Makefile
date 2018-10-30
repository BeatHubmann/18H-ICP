CXX = g++-8
CXXFLAGS = -std=c++1z -O3 -march=native -funroll-loops 

.PHONY: all
all: task figures 

task: task.cpp

results: task
	./task 32 1100000 1 0.01 > results.dat

figures: results figures.gplt
	gnuplot figures.gplt

.PHONY: clean
clean:
	@rm -rf task *.dat *.eps *.aux *.log *.synctex 