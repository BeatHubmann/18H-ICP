CXX = g++-8
CXXFLAGS = -std=c++1z -O3 -march=native -funroll-loops -ffast-math
CXXFLAGS += -I/usr/local/include/eigen3


.PHONY: all
all: main output.dat figures 

main: main.cpp 

output.dat: main
	./main 80 40.0 1.0 9.81 0.001 15000 > output.dat

.PHONY: figures
figures: output.dat figures.gplt
	gnuplot figures.gplt

.PHONY: clean
clean:
	@rm -rf main *.dat *.aux *.log *.synctex 