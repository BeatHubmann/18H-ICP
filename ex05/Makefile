CXX = g++-8
CXXFLAGS = -std=c++1z -O3 -march=native -funroll-loops 

.PHONY: all
all: task figures 

task: task.cpp

figures: task figures.gplt
	./task 1024 1 > results.dat
	gnuplot figures.gplt

.PHONY: clean
clean:
	@rm -rf task *.dat *.eps *.aux *.log *.synctex 