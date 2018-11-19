CXX = g++-8
CXXFLAGS = -std=c++1z -O3 -march=native -funroll-loops -ffast-math

.PHONY: all
all: main task1.dat task2.dat figures 

main: main.cpp 

task1.dat: main
	./main 8192 10000 1 1.0 0 0 1.0 > task1.dat

task2.dat: main
	./main 8192 10000 1 1.0 1 1 1.0 > task2.dat

.PHONY: figures
figures: task1.dat task2.dat figures.gplt
	gnuplot figures.gplt

.PHONY: clean
clean:
	@rm -rf main *.dat *.aux *.log *.synctex 